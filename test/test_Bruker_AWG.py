from autodeer.hardware import BrukerAWG

config_file = "test/test_data/test_Bruker_config.yaml"

def test_create_interface():
    interface = BrukerAWG(config_file)
