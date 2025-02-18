import csv
import glob
import gzip
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime

# Configurar logger
logger = logging.getLogger("condiga 0.2.2")
logging.basicConfig(level=logging.INFO)

def create_directory(path, clear=False):
    if clear and os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)

def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Erro ao executar comando: {command}\n{e}")

def download_genomes(taxid_list, assembly_summary, output):
    create_directory(f"{output}/Assemblies", clear=True)
    
    taxid_data = {taxid: {"url": "", "file_path": "", "present": False} for taxid in taxid_list}
    
    with open(assembly_summary) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter="\t")
        
        for row in csv_reader:
            if row and not row[0].startswith("#") and row[5] in taxid_list:
                taxid, version_status, assembly_level, genome_rep, rel_date, url = row[5], row[10], row[11], row[13], row[14].split("/"), row[19]
                
                if version_status == "latest" and genome_rep == "Full" and assembly_level in ["Complete Genome", "Contig", "Chromosome"]:
                    myurl = f"{url}/{url.split('/')[-1]}_genomic.fna.gz"
                    local_file = os.path.join(output, "Assemblies", f"{url.split('/')[-1]}_genomic.fna.gz")
                    myfile_name = local_file.replace(".gz", "")
                    
                    command = ""
                    if myurl.startswith("https:"):
                        command = f"rsync -P --copy-links --times --verbose {myurl.replace('https:', 'rsync:')} {output}/Assemblies/"
                    elif myurl.startswith("ftp:"):
                        command = f"rsync -P --copy-links --times --verbose {myurl.replace('ftp:', 'rsync:')} {output}/Assemblies/"
                    
                    if command:
                        run_command(command)
                    
                    if os.path.exists(local_file):
                        try:
                            with gzip.open(local_file, "rb") as f_in, open(myfile_name, "wb") as f_out:
                                shutil.copyfileobj(f_in, f_out)
                            os.remove(local_file)
                            taxid_data[taxid].update({"url": url, "file_path": myfile_name, "present": True})
                        except Exception as e:
                            logger.error(f"Erro ao descompactar {local_file}: {e}")
                            os.remove(local_file)
    
    return taxid_data

def get_ref_lengths(taxid_data):
    taxid_file_len = {}
    
    for taxid, info in taxid_data.items():
        if info["present"]:
            command = f"grep -v '>' {info['file_path']} | wc | awk '{{print $3-$1}}'"
            try:
                output = subprocess.check_output(command, shell=True).decode("utf-8").strip()
                taxid_file_len[taxid] = int(output)
            except Exception as e:
                logger.error(f"Erro ao obter comprimento de referência para {taxid}: {e}")
    
    return taxid_file_len

def rename_and_copy_genomes(taxid_data, species_data, output, rel_abundance, genome_coverage):
    create_directory(f"{output}/Reference_Sequences", clear=True)
    
    n_species, n_taxid = 0, 0
    
    with open(f"{output}/species_stats.tsv", "w") as myfile:
        myfile.write("Species name\tRelative abundance\tGenome coverage\n")
        
        for species, values in species_data.items():
            if values["rel_abundance"] > rel_abundance and values["genome_coverage"] > genome_coverage:
                myfile.write(f"{species}\t{values['rel_abundance']}\t{values['genome_coverage']}\n")
                n_species += 1
                
                for taxid in values["taxids"]:
                    if taxid in taxid_data and taxid_data[taxid]["present"]:
                        n_taxid += 1
                        shutil.copy(taxid_data[taxid]["file_path"], f"{output}/Reference_Sequences/{taxid}.fna")
    
    logger.info(f"{n_species} espécies identificadas com {n_taxid} taxids")

def get_ref_ids(output):
    reference_files = glob.glob(f"{output}/Reference_Sequences/*.fna")
    ref_ids = {}
    
    for ref in reference_files:
        ref_name = os.path.basename(ref).replace(".fna", "")
        command = f"grep '^>' {ref}"
        
        try:
            entries = subprocess.check_output(command, shell=True).decode("utf-8").strip().split("\n")
            ref_ids[ref_name] = [entry.split(" ")[0][1:] for entry in entries]
        except Exception as e:
            logger.error(f"Erro ao obter IDs de referência para {ref}: {e}")
    
    return ref_ids
