from extra.useful import fasta_extractor
from modules.blaster import blastn_dic, blastn_blaster, repetitive_blaster
from modules.files_manager import folder_creator

repetitive_blaster("../genome_data/ZTest2_BEDOPS_laptop/TriTrypDB-67_LinfantumJPCM5_Genome.fasta",
                   "../genome_data/ZTest2_BEDOPS_laptop/First_Blaster.csv",
                   "../genome_data/ZTest2_BEDOPS_laptop/Results/",
                   "LinJ",
                   600,
                   1,
                   20                  
)