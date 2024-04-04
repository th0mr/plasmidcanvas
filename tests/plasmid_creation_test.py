import unittest
import os
from pBR322_basic_example import get_basic_plasmid
from pBR322_overlapping_example import get_overlapping_plasmid
from pBR322_curved_example import get_curved_plasmid

class TestPlasmidCreation(unittest.TestCase):
    """
    Runs the three examples from the usage guide and puts them into a folder so the developer
    can inspect them manually.
    """

    test_dir = os.path.dirname(os.path.abspath(__file__))
    
    output_folder = os.path.join(test_dir, "test_output")

    basic_output_file =  os.path.join(output_folder, "basic_plasmid.png")
    ol_output_file =  os.path.join(output_folder, "overlapping_plasmid.png")
    curved_output_file = os.path.join(output_folder, "curved_plasmid.png")

    output_files = [basic_output_file, ol_output_file, curved_output_file]

    # Wipe all the images from the output file and create output folder if it does not exist
    @classmethod
    def setUpClass(self):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder, exist_ok=True)
        for path in self.output_files:
            if os.path.exists(path):
                os.remove(path)
        

    def test_save_basic_plasmid_to_file(self):
        basic_plasmid = get_basic_plasmid()

        # Save to file
        basic_plasmid.save_to_file(self.basic_output_file)

        # Assertions
        self.assertTrue(os.path.exists(self.basic_output_file), "Output file does not exist")

    def test_save_overlapping_plasmid_to_file(self):
        ol_plasmid = get_overlapping_plasmid()

        # Save to file
        ol_plasmid.save_to_file(self.ol_output_file)

        # Assertions
        self.assertTrue(os.path.exists(self.ol_output_file), "Output file does not exist")
        
    def test_save_curved_plasmid_to_file(self):
        curved_plasmid = get_curved_plasmid()

        # Save to file
        curved_plasmid.save_to_file(self.curved_output_file)

        # Assertions
        self.assertTrue(os.path.exists(self.curved_output_file), "Output file does not exist")

if __name__ == '__main__':
    unittest.main()