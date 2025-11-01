#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.spatial import KDTree  # Fast GEM hole lookup


class process_tracks:
    def __init__(self, infile, outpath, nGEM, drift, min_drift_length, max_drift_length, v_drift, rotate_primary_track,
                 sigmaT, sigmaL, sigmaT_trans, sigmaL_trans, sigmaT_induc, sigmaL_induc, drift_gap_length,
                 GEM_width, GEM_height, GEM_thickness, hole_diameter, hole_pitch, extra_GEM_diffusion, amplify, gain,
                 transfer_gap_length, induction_gap_length, GEM_offsetsx, GEM_offsetsy, force_through_hole,
                 randomize_position, write_ITO, sigmaTe=None, sigmaLe=None, cam_bins_x=2048,
                 cam_bins_y=1152, cam_width=8, cam_height=4.5, write_gain=False, overwrite=False, use_gpu=False):

        self.gpu = use_gpu
        self.nGEM = int(nGEM)
        self.gain = gain ** (1 / self.nGEM)  # scale of random exponential for amplification per GEM
        self.drift_gap_length = drift_gap_length

        self.cam_bins_x = cam_bins_x
        self.cam_bins_y = cam_bins_y
        self.cam_width = cam_width
        self.cam_height = cam_height

        self.write_ITO = write_ITO

        self.randomize_position = randomize_position
        self.rotate = rotate_primary_track
        self.v_drift = v_drift
        self.sigmaT = float(sigmaT)
        self.sigmaL = float(sigmaL)
        self.sigmaT_trans = float(sigmaT_trans)
        self.sigmaL_trans = float(sigmaL_trans)
        self.sigmaT_induc = float(sigmaT_induc)
        self.sigmaL_induc = float(sigmaL_induc)

        self.GEM_width = GEM_width
        self.GEM_height = GEM_height
        self.GEM_thickness = GEM_thickness
        self.hole_diameter = hole_diameter
        self.hole_pitch = hole_pitch
        self.GEM_offsetsx = GEM_offsetsx
        self.GEM_offsetsy = GEM_offsetsy
        self.extra_GEM_diffusion = extra_GEM_diffusion
        self.force_through_hole = force_through_hole
        self.sigmaTe = None if sigmaTe is None else float(sigmaTe)
        self.sigmaLe = None if sigmaLe is None else float(sigmaLe)

        print(infile)

        self.data = pd.read_feather(infile)
        self.migdal = 'ID' in self.data.columns

        # Isotropize track direction if rotate is set to True
        if self.rotate:
            self.data = self.rotate_tracks(self.data)

        # Randomize xy position if specified
        if self.randomize_position:
            self.data['xshift'] = np.random.uniform(-self.cam_width / 2, self.cam_width / 2, len(self.data))
            self.data['yshift'] = np.random.uniform(-self.cam_height / 2, self.cam_height / 2, len(self.data))
        else:
            self.data['xshift'] = 0
            self.data['yshift'] = 0

        self.data['x'] = self.data['x'] + self.data['xshift']
        self.data['y'] = self.data['y'] + self.data['yshift']

        # Fiducialize in z
        self.data['z'] = self.data['z'].apply(lambda x: x - x.mean())
        self.data['drift_length'] = np.random.uniform(min_drift_length, max_drift_length, len(self.data))  # cm
        self.data['z'] = self.data['drift_length'] + self.data['z']

        # Fiducial clipping flags and indices
        self.data['GEM_clipped'] = self.data['z'].apply(lambda x: 1 if len(np.where(x < 0)[0]) > 0 else 0)
        self.data['Cathode_clipped'] = self.data['z'].apply(lambda x: 1 if len(np.where(x > self.drift_gap_length)[0]) > 0 else 0)
        self.data['fiducial_indices'] = self.data['z'].apply(lambda x: np.where((x >= 0) & (x <= self.drift_gap_length))[0])

        # Remove nonfiducial coordinates (x,y,z)
        for col in ['x', 'y', 'z']:
            self.data[col] = self.data.apply(lambda row: row[col][row['fiducial_indices']], axis=1)

        # Keep ID in-sync (prevents ragged arrays) — only if Migdal
        if self.migdal:
            self.data['ID'] = self.data.apply(lambda row: row['ID'][row['fiducial_indices']], axis=1)

        # Apply drift+diffusion if specified
        if drift:
            xdiff, ydiff, zdiff = [], [], []
            for i in range(len(self.data)):
                track = self.data.iloc[i]
                if self.gpu:
                    xd, yd, zd = self.apply_diffusion_gpu(track['x'], track['y'], track['z'], diffusion_length=track['z'])
                else:
                    xd, yd, zd = self.apply_diffusion(track['x'], track['y'], track['z'], diffusion_length=track['z'])
                xdiff.append(xd); ydiff.append(yd); zdiff.append(zd)
            self.data['xdiff'] = xdiff
            self.data['ydiff'] = ydiff
            self.data['zdiff'] = zdiff
        else:
            self.data['xdiff'] = self.data['x']
            self.data['ydiff'] = self.data['y']
            self.data['zdiff'] = self.data['z']

        # Amplify if specified
        if amplify:
            # Populate lists of hole positions for each GEM
            if isinstance(self.GEM_offsetsx, list) or isinstance(self.GEM_offsetsy, list):
                if len(self.GEM_offsetsx) != len(self.GEM_offsetsy):
                    raise ValueError("If GEM_offsetsx or GEM_offsetsy are declared as a list, both must have the same length.")
                elif self.nGEM > len(self.GEM_offsetsx):
                    raise ValueError("Number of GEMs must be <= number of entries in GEM_offset{x,y}.")
                else:
                    self.hole_positions = [self.create_GEM_holes(offsetx=self.GEM_offsetsx[i], offsety=self.GEM_offsetsy[i]) for i in range(self.nGEM)]
            else:  # Make GEMs aligned if offsets are not a list
                self.hole_positions = [self.create_GEM_holes(offsetx=0, offsety=0) for _ in range(self.nGEM)]

            # KDTree per GEM
            self.hole_trees = [KDTree(self.hole_positions[i]) for i in range(self.nGEM)]

            # Prebuild CUDA hole positions if snapping on GPU
            if self.gpu and self.force_through_hole:
                import torch  # local import to avoid hard dependency when gpu=False
                self.hole_positions_torch = [
                    torch.tensor(np.asarray(self.hole_positions[i], dtype=np.float32), device="cuda", dtype=torch.float32)
                    for i in range(self.nGEM)
                ]

            # Output branches
            xcam, ycam, qcam, fraccam = [], [], [], []
            xITO, zITO, qITO, fracITO = [], [], [], []

            # Apply gain and digitize on a track-by-track basis
            for evt_idx in tqdm(range(len(self.data))):
                event = self.data.iloc[evt_idx]
                track = self._build_track_from_event(event)  # robust builder (no ragged arrays)

                # Loop through GEMs
                for gem_idx in range(self.nGEM):
                    if len(track) == 0:
                        break

                    # Filter or snap to hole centers
                    if not self.force_through_hole:
                        # Vectorized hole acceptance
                        dists, _ = self.hole_trees[gem_idx].query(track[:, :2])
                        mask = dists <= self.hole_radius_cm
                        track = track[mask]
                    else:
                        if self.gpu:
                            track = self.snap_to_nearest_holes_gpu(track, self.hole_positions_torch[gem_idx])
                        else:
                            track = self.snap_to_nearest_holes(track, self.hole_trees[gem_idx], self.hole_positions[gem_idx], keep_radius_gate=False)

                    if len(track) == 0:
                        break

                    gap_length = transfer_gap_length + self.extra_GEM_diffusion if gem_idx != self.nGEM - 1 else self.extra_GEM_diffusion
                    track = self.apply_amplification(track, gap_length)

                # Camera digitization
                if len(track) > 0:
                    if self.migdal:
                        xc, yc, qc, fracc = self.digitize_camera_migdal(track[:, 0], track[:, 1], track[:, 3])
                        fraccam.append(fracc)
                    else:
                        xc, yc, qc = self.digitize_camera(track[:, 0], track[:, 1])
                else:
                    if self.migdal:
                        xc, yc, qc, fracc = [], [], [], []
                        fraccam.append(fracc)
                    else:
                        xc, yc, qc = [], [], []
                xcam.append(xc)
                ycam.append(yc)
                qcam.append(qc)

                # Induction smearing + ITO digitization
                try:
                    track = self.smear_induction(track, induction_gap_length)
                    if len(track) > 0:
                        if self.migdal:
                            xi, zi, qi, fraci = self.digitize_ITO_migdal(track[:, 0], track[:, 2], track[:, 3])
                            fracITO.append(fraci)
                        else:
                            xi, zi, qi = self.digitize_ITO(track[:, 0], track[:, 2])
                    else:
                        if self.migdal:
                            xi, zi, qi, fraci = [], [], [], []
                            fracITO.append(fraci)
                        else:
                            xi, zi, qi = [], [], []
                except Exception:
                    if self.migdal:
                        xi, zi, qi, fraci = [], [], [], []
                        fracITO.append(fraci)
                    else:
                        xi, zi, qi = [], [], []
                xITO.append(xi)
                zITO.append(zi)
                qITO.append(qi)

                if self.gpu:
                    try:
                        import torch
                        torch.cuda.empty_cache()
                    except Exception:
                        pass

            # Save the digitized outputs to self.data.
            self.data['xcam'] = xcam
            self.data['ycam'] = ycam
            self.data['qcam'] = [np.array(q).astype('uint16') for q in qcam]

            if self.write_ITO:
                self.data['xITO'] = xITO
                self.data['zITO'] = zITO
                self.data['qITO'] = qITO
            if self.migdal:
                self.data['ER_frac_cam'] = fraccam
                if self.write_ITO:
                    self.data['ER_frac_ITO'] = fracITO

        # Save output dataframe
        if overwrite:
            self.data.to_feather(infile)
        else:
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            outname = os.path.split(os.path.splitext(infile)[0])[1] + '_%sGEMs_%sgain_%ssigmaT_%ssigmaTransder_digitized.feather' % (
                self.nGEM, gain, self.sigmaT, sigmaT_trans)
            self.data.to_feather(os.path.join(outpath, outname))

    # ---------- helpers / utilities ----------

    def _build_track_from_event(self, event):
        """
        Build (N,3)/(N,4) ndarray, using float32 for speed and int8 for ID.
        Prevents ragged arrays by trimming to the common length.
        """
        x = np.asarray(event['xdiff'], dtype=np.float32)
        y = np.asarray(event['ydiff'], dtype=np.float32)
        z = np.asarray(event['zdiff'], dtype=np.float32)
        n = min(len(x), len(y), len(z))

        if self.migdal:
            ID = np.asarray(event['ID'], dtype=np.int8)
            n = min(n, len(ID))
            if n == 0:
                return np.empty((0, 4), dtype=np.float32)
            xyz = np.column_stack((x[:n], y[:n], z[:n]))  # (n,3) float32
            return np.column_stack((xyz, ID[:n]))         # (n,4) last col int8
        else:
            if n == 0:
                return np.empty((0, 3), dtype=np.float32)
            return np.column_stack((x[:n], y[:n], z[:n]))  # (n,3)

    def rotate_track(self, track, init_dir):
        def random_theta_phi():
            ctheta = np.random.uniform(-1, 1)
            phi = np.random.uniform(0, 2 * np.pi)
            theta = np.arccos(ctheta)
            x = np.sin(theta) * np.cos(phi)
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(theta)
            return theta, np.arctan2(y, x)

        def rotate_y(x, y, z, angle):
            xp = np.cos(angle) * x + np.sin(angle) * z
            yp = y
            zp = -np.sin(angle) * x + np.cos(angle) * z
            return xp, yp, zp

        def rotate_z(x, y, z, angle):
            xp = np.cos(angle) * x - np.sin(angle) * y
            yp = np.sin(angle) * x + np.cos(angle) * y
            zp = z
            return xp, yp, zp

        theta, phi = random_theta_phi()

        x_r1, y_r1, z_r1 = rotate_y(track['x'], track['y'], track['z'], -(np.pi / 2 - theta))
        x_r2, y_r2, z_r2 = rotate_z(x_r1, y_r1, z_r1, phi)

        dir_rx1, dir_ry1, dir_rz1 = rotate_y(init_dir[0], init_dir[1], init_dir[2], -(np.pi / 2 - theta))
        dir_rx2, dir_ry2, dir_rz2 = rotate_z(dir_rx1, dir_ry1, dir_rz1, phi)

        dir_r = np.array([dir_rx2, dir_ry2, dir_rz2])

        return x_r2, y_r2, z_r2, dir_r

    def rotate_tracks(self, df):
        xrot, yrot, zrot, dirrot = [], [], [], []
        for i in range(len(df)):
            track = df.iloc[i]
            xr, yr, zr, dirr = self.rotate_track(track, init_dir=track['truth_dir'])
            xrot.append(xr); yrot.append(yr); zrot.append(zr); dirrot.append(dirr)
        df['x'] = xrot; df['y'] = yrot; df['z'] = zrot; df['truth_dir'] = dirrot
        return df

    # --- diffusion (CPU/GPU) ---

    def apply_diffusion(self, xs, ys, zs, diffusion_length):
        xs = np.asarray(xs, dtype=np.float32)
        ys = np.asarray(ys, dtype=np.float32)
        zs = np.asarray(zs, dtype=np.float32)
        L = np.asarray(diffusion_length, dtype=np.float32)
        x_diff = np.sqrt(L) * self.sigmaT * 1e-4 * np.random.normal(0, 1, len(zs)).astype(np.float32)
        y_diff = np.sqrt(L) * self.sigmaT * 1e-4 * np.random.normal(0, 1, len(zs)).astype(np.float32)
        z_diff = np.sqrt(L) * self.sigmaL * 1e-4 * np.random.normal(0, 1, len(zs)).astype(np.float32)
        return xs + x_diff, ys + y_diff, zs + z_diff

    def apply_diffusion_gpu(self, xs, ys, zs, diffusion_length):
        import torch
        x = torch.as_tensor(xs, device="cuda", dtype=torch.float32)
        y = torch.as_tensor(ys, device="cuda", dtype=torch.float32)
        z = torch.as_tensor(zs, device="cuda", dtype=torch.float32)
        L = torch.as_tensor(diffusion_length, device="cuda", dtype=torch.float32)

        n = z.shape[0]
        if n == 0:
            return xs, ys, zs

        x = x + torch.sqrt(L) * self.sigmaT * 1e-4 * torch.randn(n, device="cuda", dtype=torch.float32)
        y = y + torch.sqrt(L) * self.sigmaT * 1e-4 * torch.randn(n, device="cuda", dtype=torch.float32)
        z = z + torch.sqrt(L) * self.sigmaL * 1e-4 * torch.randn(n, device="cuda", dtype=torch.float32)
        return x.cpu().numpy(), y.cpu().numpy(), z.cpu().numpy()

    # --- GEM geometry ---

    def create_GEM_holes(self, offsetx=0, offsety=0):
        # Convert micrometers to centimeters
        hole_diameter_cm = self.hole_diameter / 10000.0  # cm
        self.hole_radius_cm = hole_diameter_cm / 2.0
        hole_pitch_cm = self.hole_pitch / 10000.0  # cm

        # Calculate the number of holes in the x and y directions
        num_holes_x = int(self.GEM_width / hole_pitch_cm)
        num_holes_y = int(self.GEM_height / (hole_pitch_cm * np.sqrt(3) / 2))

        # Generate positions
        hole_positions = []
        for i in range(num_holes_x):
            for j in range(num_holes_y):
                x = i * hole_pitch_cm
                y = j * (hole_pitch_cm * np.sqrt(3) / 2)
                if j % 2 == 1:
                    x += hole_pitch_cm / 2
                x += offsetx / 10000.0
                y += offsety / 10000.0
                hole_positions.append((x, y))

        shift = self.GEM_width / 2.0
        pos = pd.Series(hole_positions).apply(np.array) - shift
        hole_positions = pos[pos.apply(lambda v: (v[0] > -shift) & (v[1] > -shift) & (v[0] < shift) & (v[1] < shift))].apply(tuple).to_list()
        return hole_positions

    # --- hole snapping (CPU & GPU) ---

    def snap_to_nearest_holes(self, track, tree, hole_positions, keep_radius_gate=False):
        if track is None or len(track) == 0:
            return track
        xy = track[:, :2]
        dist, idx = tree.query(xy)
        if keep_radius_gate:
            mask = dist <= self.hole_radius_cm
            return track[mask]
        track[:, :2] = np.asarray(hole_positions, dtype=np.float32)[idx]
        return track

    def snap_to_nearest_holes_gpu(self, track, hole_positions_t):
        if track is None or len(track) == 0:
            return track
        import torch
        t = torch.as_tensor(track, device="cuda", dtype=torch.float32)
        xy = t[:, :2]                               # (N,2)
        d = torch.cdist(xy, hole_positions_t)       # (N, M)
        idx = torch.argmin(d, dim=1)                # (N,)
        snapped_xy = hole_positions_t[idx]          # (N,2)
        t[:, :2] = snapped_xy
        return t.cpu().numpy()

    # --- gain / amplification ---

    def generate_gain_points(self, x, x_post, gain_electrons, gap_length, diff_coeff, extra_sigma):
        # CPU: add diffusion (gap) + optional extra_sigma (e.g., σTe/σLe) as independent Gaussians
        for enum, val in np.ndenumerate(gain_electrons):
            start_ind = int(np.sum(gain_electrons[:enum[0]]))
            end_ind = int(np.sum(gain_electrons[:enum[0] + 1]))
            s_gap = np.sqrt(gap_length) * diff_coeff * 1e-4
            if extra_sigma is not None:
                s_extra = extra_sigma * 1e-4
                sigma_total = np.sqrt(s_gap**2 + s_extra**2)
            else:
                sigma_total = s_gap
            x_post[start_ind:end_ind] = x[enum] + sigma_total * np.random.normal(0, 1, val).astype(np.float32)

    def generate_gain_points_GPU(self, x, gain_electrons, sigma_total):
        """
        x: 1D numpy (float32), gain_electrons: torch.int64 counts per primary,
        sigma_total: float (float32) total sigma for this axis (gap + extra) in cm.
        """
        import torch
        ge = gain_electrons.to(device='cuda')
        xv = torch.as_tensor(x, device='cuda', dtype=torch.float32)
        total = int(torch.sum(ge).item())
        if total == 0:
            return torch.empty((0,), device='cuda', dtype=torch.float32)
        noise = sigma_total * torch.randn(total, device='cuda', dtype=torch.float32)
        x_rep = xv.repeat_interleave(ge)
        return x_rep + noise

    def GEM_gain_and_diffusion(self, x, y, z, gap_length, diff_coeff_trans, diff_coeff_long):
        if not self.gpu:
            gain_electrons = np.random.exponential(self.gain, len(x)).astype(np.int64)
            x_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)), dtype=np.float32)
            y_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)), dtype=np.float32)
            z_post = np.ascontiguousarray(np.zeros(np.sum(gain_electrons)), dtype=np.float32)

            self.generate_gain_points(x, x_post, gain_electrons, gap_length=gap_length, diff_coeff=diff_coeff_trans, extra_sigma=self.sigmaTe)
            self.generate_gain_points(y, y_post, gain_electrons, gap_length=gap_length, diff_coeff=diff_coeff_trans, extra_sigma=self.sigmaTe)
            self.generate_gain_points(z, z_post, gain_electrons, gap_length=gap_length, diff_coeff=diff_coeff_long,  extra_sigma=self.sigmaLe)

            return x_post, y_post, z_post
        else:
            import torch
            # match numpy exponential with same scale (self.gain)
            exponential_dist = torch.distributions.Exponential(1.0 / self.gain)
            gain_electrons = exponential_dist.sample((len(x),)).to(dtype=torch.int64, device='cuda')

            # Base sigma from gap diffusion (cm)
            sT_gap = (gap_length ** 0.5) * diff_coeff_trans * 1e-4
            sL_gap = (gap_length ** 0.5) * diff_coeff_long  * 1e-4
            # Extra blur terms (cm) if provided
            sTe = (self.sigmaTe * 1e-4) if (self.sigmaTe is not None) else 0.0
            sLe = (self.sigmaLe * 1e-4) if (self.sigmaLe is not None) else 0.0
            # Combine in quadrature (independent Gaussians)
            sX = np.float32(np.sqrt(sT_gap**2 + sTe**2))
            sY = np.float32(np.sqrt(sT_gap**2 + sTe**2))
            sZ = np.float32(np.sqrt(sL_gap**2 + sLe**2))

            x_post = self.generate_gain_points_GPU(x, gain_electrons, sX)
            y_post = self.generate_gain_points_GPU(y, gain_electrons, sY)
            z_post = self.generate_gain_points_GPU(z, gain_electrons, sZ)

            return x_post.cpu().numpy(), y_post.cpu().numpy(), z_post.cpu().numpy()

    def apply_amplification(self, track, gap_length):
        if self.migdal:
            # Separate out NR and ER charges.
            NRidx = np.where(track[:, 3] == 1)[0]
            ERidx = np.where(track[:, 3] == 0)[0]

            if len(NRidx) > 0:
                xNR, yNR, zNR = self.GEM_gain_and_diffusion(
                    track[NRidx, 0], track[NRidx, 1], track[NRidx, 2],
                    gap_length=gap_length,
                    diff_coeff_trans=self.sigmaT_trans,
                    diff_coeff_long=self.sigmaL_trans
                )
            else:
                xNR, yNR, zNR = np.array([], dtype=np.float32), np.array([], dtype=np.float32), np.array([], dtype=np.float32)

            if len(ERidx) > 0:
                xER, yER, zER = self.GEM_gain_and_diffusion(
                    track[ERidx, 0], track[ERidx, 1], track[ERidx, 2],
                    gap_length=gap_length,
                    diff_coeff_trans=self.sigmaT_trans,
                    diff_coeff_long=self.sigmaL_trans
                )
            else:
                xER, yER, zER = np.array([], dtype=np.float32), np.array([], dtype=np.float32), np.array([], dtype=np.float32)

            # Recombine & rebuild track with IDs (NR=1, ER=0)
            new_x = np.concatenate((xNR, xER))
            new_y = np.concatenate((yNR, yER))
            new_z = np.concatenate((zNR, zER))
            new_ID = np.concatenate((np.ones(len(xNR), dtype=np.int8), np.zeros(len(xER), dtype=np.int8)))
            out = np.column_stack((new_x, new_y, new_z, new_ID))
            return out
        else:
            xgain, ygain, zgain = self.GEM_gain_and_diffusion(
                track[:, 0], track[:, 1], track[:, 2],
                gap_length=gap_length,
                diff_coeff_trans=self.sigmaT_trans,
                diff_coeff_long=self.sigmaL_trans
            )
            return np.column_stack((xgain, ygain, zgain))

    # --- induction smearing (CPU/GPU) ---

    def smear_induction(self, track, induction_gap_length):
        if track is None or len(track) == 0:
            return track
        sT = (induction_gap_length ** 0.5) * self.sigmaT_induc * 1e-4
        sL = (induction_gap_length ** 0.5) * self.sigmaL_induc * 1e-4
        if not self.gpu:
            track[:, 0] += sT * np.random.normal(0, 1, len(track)).astype(np.float32)
            track[:, 1] += sT * np.random.normal(0, 1, len(track)).astype(np.float32)
            track[:, 2] += sL * np.random.normal(0, 1, len(track)).astype(np.float32)
            return track
        else:
            import torch
            t = torch.as_tensor(track, device="cuda", dtype=torch.float32)
            n = t.shape[0]
            if n == 0:
                return track
            t[:, 0] += sT * torch.randn(n, device="cuda", dtype=torch.float32)
            t[:, 1] += sT * torch.randn(n, device="cuda", dtype=torch.float32)
            t[:, 2] += sL * torch.randn(n, device="cuda", dtype=torch.float32)
            return t.cpu().numpy()

    # --- digitizers ---

    def digitize_camera_migdal(self, x, y, ID):
        NRidx = np.where(ID == 1)[0]
        ERidx = np.where(ID == 0)[0]
        NRhist = np.histogram2d(x[NRidx], y[NRidx], bins=(self.cam_bins_x, self.cam_bins_y),
                                range=((-self.cam_width / 2, self.cam_width / 2),
                                       (-self.cam_height / 2, self.cam_height / 2)))[0].T
        ERhist = np.histogram2d(x[ERidx], y[ERidx], bins=(self.cam_bins_x, self.cam_bins_y),
                                range=((-self.cam_width / 2, self.cam_width / 2),
                                       (-self.cam_height / 2, self.cam_height / 2)))[0].T
        totalhist = NRhist + ERhist
        fraction = np.divide(ERhist, totalhist, out=np.zeros_like(ERhist, dtype=float), where=totalhist != 0)
        sparse_hist = np.where(totalhist > 0)
        y_idx, x_idx = sparse_hist
        q = totalhist[sparse_hist]
        frac = fraction[sparse_hist]
        return x_idx.astype(np.int32), y_idx.astype(np.int32), q.astype(np.uint16), frac.astype(np.float32)

    def digitize_camera(self, x, y):
        a = np.histogram2d(x, y, bins=(self.cam_bins_x, self.cam_bins_y),
                           range=((-self.cam_width / 2, self.cam_width / 2),
                                  (-self.cam_height / 2, self.cam_height / 2)))[0].T
        sparse_hist = np.where(a > 0)
        y_idx, x_idx = sparse_hist
        q = a[sparse_hist]
        return x_idx.astype(np.int32), y_idx.astype(np.int32), q.astype(np.uint16)

    def digitize_ITO_migdal(self, x, z, ID):
        NRidx = np.where(ID == 1)[0]
        ERidx = np.where(ID == 0)[0]
        NRhist = np.histogram2d(x[NRidx], z[NRidx], bins=(120, 150), range=((-5, 5), (0, 3.9)))[0].T
        ERhist = np.histogram2d(x[ERidx], z[ERidx], bins=(120, 150), range=((-5, 5), (0, 3.9)))[0].T
        totalhist = NRhist + ERhist
        fraction = np.divide(ERhist, totalhist, out=np.zeros_like(ERhist, dtype=float), where=totalhist != 0)
        sparse_hist = np.where(totalhist > 0)
        z_idx, x_idx = sparse_hist
        q = totalhist[sparse_hist]
        frac = fraction[sparse_hist]
        return x_idx.astype(np.int32), z_idx.astype(np.int32), q.astype(np.uint16), frac.astype(np.float32)

    def digitize_ITO(self, x, z):
        a = np.histogram2d(x, z, bins=(120, 150), range=((-5, 5), (0, 3.9)))[0].T
        sparse_hist = np.where(a > 0)
        z_idx, x_idx = sparse_hist
        q = a[sparse_hist]
        return x_idx.astype(np.int32), z_idx.astype(np.int32), q.astype(np.uint16)


# === Entry point ===
if __name__ == '__main__':
    import yaml
    from os import sys

    """Load configuration.yaml"""
    with open('configuration.yaml', 'r') as cfg:
        config = yaml.safe_load(cfg)
        settings = config['Sim_settings']
        gas_cfg = config['Gas_props']
        tpc_cfg = config['TPC_sim']

    # Handle command line override
    infile = sys.argv[1] if len(sys.argv) > 1 else settings['digitization_input_file']

    # Unpack config values
    drift = settings['apply_drift']
    rotate = settings['rotate_tracks']
    nGEM = tpc_cfg['nGEM']
    min_drift_length = tpc_cfg['min_drift_length']
    max_drift_length = tpc_cfg['max_drift_length']
    v_drift = gas_cfg['vd']
    sigmaT = gas_cfg['sigmaT']
    sigmaL = gas_cfg['sigmaL']
    sigmaT_trans = gas_cfg['sigmaT_trans']
    sigmaL_trans = gas_cfg['sigmaL_trans']
    sigmaT_induc = gas_cfg['sigmaT_induc']
    sigmaL_induc = gas_cfg['sigmaL_induc']
    drift_gap_length = tpc_cfg['drift_gap_length']
    GEM_width = tpc_cfg['GEM_width']
    GEM_height = tpc_cfg['GEM_height']
    GEM_thickness = tpc_cfg['GEM_thickness']
    hole_diameter = tpc_cfg['hole_diameter']
    hole_pitch = tpc_cfg['hole_pitch']
    trans_gap = tpc_cfg['transfer_gap_length']
    induc_gap = tpc_cfg['induction_gap_length']
    GEM_offsetsx = tpc_cfg['GEM_offsetsx']
    GEM_offsetsy = tpc_cfg['GEM_offsetsy']
    extra_GEM_diffusion = tpc_cfg['extra_GEM_diffusion']
    force_through_hole = tpc_cfg['force_through_GEM_hole']

    # Derived optical resolution (μm → cm handled in code paths)
    sigmaTe = np.sqrt((hole_pitch / np.sqrt(12)) ** 2 + (hole_diameter / 4) ** 2)
    sigmaLe = None  # keep None unless you want an extra long. blur in gain

    cam_bins_x = tpc_cfg['cam_bins_x']
    cam_bins_y = tpc_cfg['cam_bins_y']
    cam_width = tpc_cfg['cam_width']
    cam_height = tpc_cfg['cam_height']

    amplify = settings['apply_amplification']
    gain = tpc_cfg['gain']
    randomize_position = settings['randomize_position']
    write_gain = settings['write_gain']
    write_ITO = settings['write_ITO']
    overwrite = settings['overwrite_output']
    outpath = settings['output_dir']
    gpu = settings['gpu']

    if gpu:
        try:
            import torch
            if not torch.cuda.is_available():
                print("WARNING: GPU requested but CUDA not available. Falling back to CPU.")
                gpu = False
        except Exception:
            print("WARNING: GPU requested but PyTorch not available. Falling back to CPU.")
            gpu = False

    process_tracks(
        infile=infile,
        outpath=outpath,
        nGEM=nGEM,
        drift=drift,
        min_drift_length=min_drift_length,
        max_drift_length=max_drift_length,
        v_drift=v_drift,
        rotate_primary_track=rotate,
        sigmaT=sigmaT,
        sigmaL=sigmaL,
        sigmaT_trans=sigmaT_trans,
        sigmaL_trans=sigmaL_trans,
        sigmaT_induc=sigmaT_induc,
        sigmaL_induc=sigmaL_induc,
        sigmaTe=sigmaTe,
        sigmaLe=sigmaLe,
        drift_gap_length=drift_gap_length,
        GEM_width=GEM_width,
        GEM_height=GEM_height,
        GEM_thickness=GEM_thickness,
        hole_diameter=hole_diameter,
        hole_pitch=hole_pitch,
        transfer_gap_length=trans_gap,
        induction_gap_length=induc_gap,
        GEM_offsetsx=GEM_offsetsx,
        GEM_offsetsy=GEM_offsetsy,
        extra_GEM_diffusion=extra_GEM_diffusion,
        force_through_hole=force_through_hole,
        cam_bins_x=cam_bins_x,
        cam_bins_y=cam_bins_y,
        cam_width=cam_width,
        cam_height=cam_height,
        randomize_position=randomize_position,
        amplify=amplify,
        write_gain=write_gain,
        write_ITO=write_ITO,
        gain=gain,
        overwrite=overwrite,
        use_gpu=gpu,
    )
