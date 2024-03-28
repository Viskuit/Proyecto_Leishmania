from extra.useful import fasta_extractor
from modules.blaster import blastn_dic, blastn_blaster, repetitive_blaster
from modules.files_manager import folder_creator
import os
# print(os.getcwd())
# os.chdir("./SIDER_RepetitiveSearcher")

repetitive_blaster("../../genome_data/ZTest3_laptop_with_chr32/TriTrypDB-67_LinfantumJPCM5_Genome.fasta",
                   "../../genome_data/ZTest3_laptop_with_chr32/First_Blaster.csv",
                   "../../genome_data/ZTest3_laptop_with_chr32/Results/",
                   "LinJ",
                   600,
                   2,
                   20                  
)