
### STEP2. Run GCTA ###

gcta=/SAN/mottlab/heterosis/bin/gcta/gcta64
plink=/home/ucbtlcu/bin/plink2
wd=/SAN/mottlab/heterosis/1_F2pigs/5VC/2_921Inds
datdir=/SAN/mottlab/heterosis/1_F2pigs/5VC/data

cd $wd

traits=(`cat $datdir/F2_resrm.traitnames`)
let N=${#traits[@]}-1 
idx=(`seq 0 $N`)

$plink --noweb --tfile $datdir/F2_res --maf 0.05 --geno 0.1 --remove $datdir/F2_Info_RemoveInd.txt --make-bed --out F2_Final

$gcta --bfile F2_Final --thread-num 10  --autosome --make-grm-gz --out F2_Final_add
$gcta --bfile F2_Final --thread-num 10  --autosome --make-grm-d-gz --out F2_Final_dom

$gcta --bfile F2_Final --thread-num 10  --autosome --make-grm --out F2_Final_add
$gcta --bfile F2_Final --thread-num 10  --autosome --make-grm-d --out F2_Final_dom

mkdir add_dom_VC

for i in ${idx[@]};do
    let j=$i+1
    $gcta --reml --mgrm $datdir/add_domi.txt --pheno $datdir/F2_resrm.pheno --mpheno $j --out ./add_dom_VC/${traits[$i]}
done

