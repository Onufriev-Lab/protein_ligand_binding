#!/bin/bash
set -x
debug="True"

### 1) generate crude truncated receptor for residue selection
#run in python 3.9, I used miniconda to create a custom environment
#Dependencies: AmberTools25
#usage: bashTruncAndCapLinux.sh <pdb code> <ligand name> <trunc radius>

#script will generate trunc.pdb and trunc_cap.pdb. trunc.pdb should be fine, trunc_cap.pdb breaks sometimes. 
LigandRes=$2
structure=$1
complex_original_pdb=${structure}.pdb

Radius=$3

remove_single_residues="False"
regex_numeric='^[0-9]+$'

# clean pdb
grep -E '^ATOM|^HETATM' "$complex_original_pdb" > complex_orig.pdb

## extract cpptraj residue list
cpptraj -p complex_orig.pdb --resmask \* > cpptraj_reslist.txt


#list of residues with no waters
cpptraj -p complex_orig.pdb --resmask '!:WAT,HOH' > cpptraj_reslistnohoh.txt


LigandResName=$LigandRes

# Check if column5 is numeric (to detect if it's a chain ID or residue number)
starting_residue=$(awk '/^(ATOM|HETATM)/ {print substr($0, 23, 4); exit}' complex_orig.pdb)
echo "$starting_residue"

newres=$(awk '/^(ATOM|HETATM)/ {residue = substr($0, 18, 3); num = substr($0, 23, 4); print residue, num}' complex_orig.pdb | sort -u -n -k2,2)
echo "$starting_residue"
echo "HELLO"
#echo "$newres" > grepped.txt

last_residue=$(awk '/^(ATOM|HETATM)/ {residue = substr($0, 18, 3); num = substr($0, 23, 4)}' complex_orig.pdb \
                | sort -u -n -k2,2 \
                | grep -E '^GLY|^ALA|^VAL|^LEU|^ILE|^THR|^SER|^MET|^CYS|^PRO|^PHE|^TYR|^TRP|^HIS|^LYS|^ARG|^ASP|^GLU|^ASN|^GLN' \
                | awk 'END {print $2}')
chainID=$(awk '/^ATOM/ {print substr($0, 22, 1)}' complex_orig.pdb | sort -u)

# generate truncated protein
cat <<eof> trunc.in
parm $complex_original_pdb
reference $complex_original_pdb
trajin $complex_original_pdb
strip !(:$LigandResName<:$Radius)
#strip :$LigandResName
strip :WAT,HOH
strip :Na+
strip :Cl-
trajout trunc.pdb
run
quit
eof

cpptraj -i trunc.in -o trunc.out

### 2) create the list of protein residues present in the truncated version
grep 'ATOM' trunc.pdb | awk '{print substr($0, 23, 4)}' | uniq | awk 'NR==1{first=$1;last=$1;next} $1 == last+1 {last=$1;next} {print first,last;first=$1;last=first} END{print first,last}' > resNo.txt

### 3) we need to patch the residue list to add extra residues which will then be mutated into Me Caps

## 3.1) clean initial residue list (concatenate sequences where only 1 residue separates the peptides, eg 183 and 185)

# calculate differences between ending and starting residues
awk 'BEGIN { OFS = FS } NR == 1 { last = $2; print;  next }{ $3 = $1 - last; last = $2  }1' resNo.txt > resNo_diff.txt

# get lines where difference between end and start < 2 
awk '{if ($3 > 0 && $3 <= 2 && $3 != "") { print 1 } else { print 0 } }' resNo_diff.txt > lines_to_patch.txt

# patch residue list so gaps smaller than 2 residues are bridged by the deleted residue
# this should also work with subsequent gaps,  eg 183 185 and 187 in the example 
if [ -f resNo_final.txt ] ; then
    rm resNo_final.txt
fi

lines=$(cat lines_to_patch.txt | wc -l )

for ((i=2 ; i<=$lines ; i+=1)); do
    lineNo=$((i-1))
    bool=$(awk -v var=$i 'NR==var { print $1} ' lines_to_patch.txt) 
    if [ $bool == 0 ] ; then
        awk -v var=$((lineNo)) 'NR==var { print $0 }' resNo.txt >> resNo_final.txt
    elif [ $bool == 1 ] ; then
        startRes_curr=$(awk -v var=$((lineNo)) 'NR==var { print $1 }' resNo.txt)
        bool_curr=1
        while [ $bool_curr == 1 ] ; do
            lineNo=$((lineNo+1))
            bool_curr=$(awk -v var=$((lineNo)) 'NR==var { print $1 }' lines_to_patch.txt)
        done
        endRes_curr=$(awk -v var=$((lineNo-1)) 'NR==var { print $2 }' resNo.txt)
        echo "$startRes_curr $endRes_curr" >> resNo_final.txt
        i=$lineNo
    fi
done

tail -n 1 resNo.txt  >> resNo_final.txt

# use final filter to remove single residues
if [ $remove_single_residues == "True" ] ; then
    awk '!($1==$2)'  resNo_final.txt > resNo_final_tmp.txt && mv resNo_final_tmp.txt resNo_final.txt
fi

## 3.2) create STRIP string for residues within distance threshold
while IFS=" " read -r startRes endRes ; do
    if [ $startRes == $endRes ] ; then
        startRes=$(awk -v ligandres_awk=$startRes '$6==ligandres_awk {print $1; exit}' cpptraj_reslist.txt)
        Res="$startRes,$Res"
    else
        startRes=$(awk -v ligandres_awk=$startRes '$6==ligandres_awk {print $1; exit}' cpptraj_reslist.txt)
        endRes=$(awk -v ligandres_awk=$endRes '$6==ligandres_awk {print $1; exit}' cpptraj_reslist.txt)
        Res="$startRes-$endRes,$Res"
    fi
done < <(tac resNo_final.txt)

Res=$(echo $Res | sed 's/,$//')

## 3.3) create STRIP string for Nterm residues
startRes=$(awk '{print $1}' resNo_final.txt | tac)
for res in $startRes ; do
    NtermRes=$((res-1))
    NtermRes=$(awk -v ligandres_awk=$NtermRes '$6==ligandres_awk {print $1; exit}' cpptraj_reslist.txt)
    AceRes="$NtermRes,$AceRes"
done
AceRes=$(echo $AceRes | sed 's/,$//')

## 3.4) create STRIP string for Cterm residues
endRes=$(awk '{print $2}' resNo_final.txt | tac)
for res in $endRes ; do
    CtermRes=$((res+1))
    CtermRes=$(awk -v ligandres_awk=$CtermRes '$6==ligandres_awk {print $1; exit}' cpptraj_reslist.txt)
    NmeRes="$CtermRes,$NmeRes"
done
NmeRes=$(echo $NmeRes | sed 's/,$//')

## 3.5) combine the necessary STRIP strings for a cpptraj input file and generate extended truncated pdb
cat <<eof> trunc_cap.in
parm $complex_original_pdb
reference $complex_original_pdb
trajin $complex_original_pdb
strip !(:$AceRes@C,O|:$NmeRes@N,H|:$Res|:$LigandResName)
strip :WAT,HOH
strip :Na+
strip :Cl-
trajout trunc_cap.pdb
run
quit
eof

cpptraj -i trunc_cap.in -o trunc_cap.out

## 3.6) rename Nterm residues to ACE
for res in $startRes ; do
    NtermRes=$((res-1))
    oldResName=$(awk -v var="$NtermRes" '{ 
        temp = substr($0, 23, 4); 
        gsub(/ /, "", temp); 
        if (temp == var) 
            print substr($0, 18, 4) 
    }' trunc_cap.pdb | tail -n 1 | tr -d ' ')
    # Updated for Linux: remove the mac-specific empty backup extension parameter in sed
    sed -E -i "s/$oldResName([[:space:]]+)(($chainID[[:space:]]+)?$NtermRes)/ACE\1\2/g" trunc_cap.pdb
done

## 3.7) rename Cterm residues to NME
for res in $endRes ; do
    CtermRes=$((res+1))
    oldResName=$(awk -v var="$CtermRes" '{ 
        temp = substr($0, 23, 4); 
        gsub(/ /, "", temp); 
        if (temp == var) 
            print substr($0, 18, 4) 
    }' trunc_cap.pdb | tail -n 1 | tr -d ' ')
    # Updated for Linux: remove the mac-specific empty backup extension parameter in sed
    sed -E -i "s/$oldResName([[:space:]]+)(($chainID[[:space:]]+)?$CtermRes)/NME\1\2/g" trunc_cap.pdb
done

# grep '^ATOM' trunc.pdb | awk '{print substr($0, 23, 4)}' | sort -u | wc -l
# echo $structure
#cat trunc.pdb > "$structure-$Radius.pdb"



#### 4) remove unnecessary files
if [ $debug == "False" ] ; then
    rm complex_orig.pdb trunc.in trunc.out resNo.txt resNo_diff.txt lines_to_patch.txt # resNo_final.txt may also be removed if desired
fi
