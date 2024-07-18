import ruamel.yaml #library preserves structure and comments of yaml file when modifying
import subprocess

def read_and_modify_yaml(yaml_file, n_tracks, energy, parallel_chunks):
    yaml = ruamel.yaml.YAML()
    
    # Read the YAML file
    with open(yaml_file, 'r') as file:
        config = yaml.load(file)

    # Modify the specified parameters
    config['Degrad_card']['n_tracks'] = n_tracks
    config['Degrad_card']['energy'] = energy
    config['Sim_settings']['parallel_chunks'] = parallel_chunks

    # Write the updated YAML content back to the file
    with open(yaml_file, 'w') as file:
        yaml.dump(config, file)

def run_processing_script(script_path):
    result = subprocess.run(['python', script_path], capture_output=True, text=True)
    if result.returncode == 0:
        print("Script executed successfully")
    else:
        print("Script execution failed")
        print(result.stdout)
        print(result.stderr)

if __name__ == "__main__":
    yaml_file_path = 'configuration.yaml'
    processing_script_path = 'run_degrad.py'

    # Modify these values as needed
    n_tracks = 400
    energies = [200 * i for i in range(10, 60)]  # 2000-12000 eV
    parallel_chunks = 1

    for e in energies:
        read_and_modify_yaml(yaml_file_path, n_tracks, e, parallel_chunks)
        run_processing_script(processing_script_path)
