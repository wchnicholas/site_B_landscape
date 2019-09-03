using CSV
using DataFrames
#names = ["HK68", "Bk79", "Bei89", "Mos99", "Bris07P194", "Bris07L194", "NDako16"]
names = ["HK68", "Bk79", "Bei89", "Mos99", "Bris07L194", "NDako16"]
#wt = ["KGSESV", "EESENV", "EEYENV", "KESDNA", "HKFDFA", "HNSDFA", "HKFDFA"]
wt = ["KGSESV", "EESENV", "EEYENV", "KESDNA", "HNSDFA", "HKFDFA"]
pos = [156, 158, 159, 190, 193, 196]
# 4 x 4 x 3 x 2 x 3 x 2 = 576
# nic said Mos99 has no wt, background four is QKYDST, but that doesn't appear either
# chose one with largest selected count
function read(n) 
    d = CSV.read("result/EpiB_Index_$n.tsv", delim="\t")
    sort!(d, :ID)
    d[:strain] = n
    d[:selected] = d[:Rep1Count] .+ d[:Rep2Count]
    d[:v] = 1./(1+d[:selected]) .+ 1./(1+d[:InputCount])
    d[:f] = log.((d[:selected]+1)./(d[:InputCount]+1))
    return d
end
#data[:f] = log.(data[:Rep1Count] .+ data[:Rep2Count])./(2*data[:InputCount])
data = map(read, names)
dall = reduce(vcat, data)
CSV.write("result/data_all.csv", dall)
