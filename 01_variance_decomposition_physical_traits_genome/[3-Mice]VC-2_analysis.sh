
### STEP2. Run GCTA ###

gcta=/SAN/mottlab/heterosis/bin/gcta/gcta64
plink=/home/ucbtlcu/bin/plink2
wd=/SAN/mottlab/heterosis/4_HSmice/5VC/AllInds
datdir=/SAN/mottlab/heterosis/4_HSmice/5VC/data

cd $wd

traits=(`cat $datdir/HSmice_res.traitnames`)
let N=${#traits[@]}-1 
idx=(`seq 0 $N`)

$plink --noweb --tfile $datdir/HSmice_res --maf 0.05 --geno 0.1 --make-bed --out HSmice_Final

$gcta --bfile HSmice_Final --thread-num 10  --autosome --make-grm-gz --out HSmice_Final_add
$gcta --bfile HSmice_Final --thread-num 10  --autosome --make-grm-d-gz --out HSmice_Final_dom

$gcta --bfile HSmice_Final --thread-num 10  --autosome --make-grm --out HSmice_Final_add
$gcta --bfile HSmice_Final --thread-num 10  --autosome --make-grm-d --out HSmice_Final_dom

mkdir add_dom_VC

for i in ${idx[@]};do
    let j=$i+1
    $gcta --reml --mgrm $datdir/add_domi.txt --pheno $datdir/HSmice_res.pheno --mpheno $j --out ./add_dom_VC/${traits[$i]}
done




