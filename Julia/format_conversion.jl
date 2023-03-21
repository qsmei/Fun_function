using Base.Threads
using LinearAlgebra
using DelimitedFiles
using StatsBase
using BenchmarkTools
using CSV,DataFrames

num_file="/usr/home/qgg/qumei/Saija/jersey_test_genotypes.txt"
map_file="/usr/home/qgg/qumei/Saija/snp_info_JER.txt"

snp_num=CSV.read(num_file,DataFrame,header=false,ntasks=Threads.nthreads());
snp_map=DelimitedFiles.readdlm(map_file,skipstart=1); #Marker chr  position  Chip1

function num_AGCT(num)  #return tuple is much faster than return vector(array),and also include the usage of memory
		 if num==0			
			return ("A","A")		
		elseif num==1   #no phasing information		 
			return ("A","T")		 
		else 		
		   return ("T","T")
		end   
end


function fimpute_to_plink(snp_num::DataFrame)

		 n_IND=size(snp_num,1)
		 n_SNP=size(snp_num,2)-3	
		 ped_AGCT=Array{String1}(undef, n_IND,n_SNP*2) #only store SNP
		 ped_id=Array{Int64}(undef, n_IND,6)	 
		 #convert snp_num to ped 
		@threads  for j in 1:n_SNP					
						for i in 1:n_IND			
							tmp=(num_AGCT(snp_num[i,j+3]))	
							ped_AGCT[i,(2*j-1)]=tmp[1];	
							ped_AGCT[i,(2*j)]=tmp[2];	
						end 
		     
				end

		 return ped_AGCT;
end

function fimpute_to_plink3(snp_num::DataFrame)

		 n_IND=size(snp_num,1)
		 n_SNP=size(snp_num,2)-3	
		 ped_AGCT=Array{String1}(undef, n_IND,n_SNP*2) #only store SNP
		 ped_id=Array{Int64}(undef, n_IND,6)	 
		 #convert snp_num to ped 
		@threads  for j in 1:n_SNP					
						for i in 1:n_IND	
							tmp=(num_AGCT1(snp_num[i,j+3]))								
							ped_AGCT[i,(2*j-1)]=tmp[1];			
							ped_AGCT[i,(2*j)]=tmp[2];	
						end 
		     
				end

		 return ped_AGCT;
end


@btime fimpute_to_plink(snp_num);
@btime fimpute_to_plink1(snp_num);
@btime fimpute_to_plink2(snp_num);
@btime fimpute_to_plink3(snp_num);


@btime for i in 1:2
	   c[1]=a[1]
	   c[2]=a[2]
	   end

@btime c[1:2]=b
   
@btime  tmp=collect(a);for i in 1:2
		c[i]=tmp[i]
		end
	 