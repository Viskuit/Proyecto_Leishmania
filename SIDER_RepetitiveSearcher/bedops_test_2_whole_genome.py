from extra.useful import fasta_extractor
from modules.blaster import blastn_dic, blastn_blaster, repetitive_blaster
from modules.files_manager import folder_creator

repetitive_blaster("../SIDER_Data/Z_BEDOPS_Test2_whole_genome/TriTrypDB-67_LinfantumJPCM5_Genome.fasta",
                   "../SIDER_Data/Z_BEDOPS_Test2_whole_genome/First_Blaster.csv",
                   "../SIDER_Data/Z_BEDOPS_Test2_whole_genome/Results/",
                   "LinJ",
                   600,
                   1,
                   2                  
)