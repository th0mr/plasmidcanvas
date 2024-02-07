

from plasmid import Plasmid
from Bio import SeqIO

def gb_file_to_plasmid(path_to_gb_file) -> Plasmid:
    gb_obj = SeqIO.read(path_to_gb_file, "genbank")

    # Extract information about the plasmid from the gb file
    plasmid_name = gb_obj.name
    plasmid_sequence = gb_obj.seq
    plasmid_length = len(plasmid_sequence)
    plasmid_annotations = gb_obj.annotations
    plasmid_features = gb_obj.features
    plasmid_feature_types=set([ feature.type for feature in plasmid_features ])

    # Debugging
    print(plasmid_name)
    print(plasmid_length)
    print(plasmid_sequence)
    print(plasmid_annotations)
    print(plasmid_features)
    print(plasmid_feature_types)

gb_file_to_plasmid("sequence.gb")