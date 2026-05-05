# Day 1 Handout

## 1. Understanding Linux Server & Installing Miniconda

• What is a Linux server  
• How to connect to a remote server from your laptop using VS Code application  
• Basic idea of working in a terminal  

• What is Miniconda and why we use it?  
• Concept of environments (like separate workspaces for different tools)  
• Create two environments for FastQC and Trimmomatic tools  

---

## 2. Installing FastQC and Trimmomatic tools for Quality Check

• Installing FastQC (quality checking tool) using conda  
• Installing Trimmomatic (read trimming tool) using conda  

---

## 3. Understanding Raw Sequence data and file formats 

• 16s Sequencing data from IBD patients  
• Number of reads and length of reads  
• Structure of FASTQ file and FASTA file  
• Conversion of file formats  

---

## 4. Quality Check of Raw Data using FastQC

• Run FastQC on all the samples  
• Understand basic interpretation of result  

---

## 5. Trimmomatic Low-Quality Reads using Trimmomatic (Paired end data)

• Run the Trimmomatic for one sample using specific parameters  
• Understand basic interpretation of result (Understand how many sequences have been preserved)  
• Run the Trimmomatic on all samples using bash script  

---

## 6. Taxonomic Classification using SPINGO for 16s Sequence data

• Install the SPINGO (Directly from GitHub)  
• Download the RDP Database  
• Run the SPINGO on one sample  
• Understand basic interpretation of result  
• Run the SPINGO on all the samples using bash script  

---

## WorkFlow
16s Raw Sequencing Data >> Quality Check >> Trimming >> Taxonomic Classification
