using Base: root_module_exists
using Distributed
using BenchmarkTools
using Dates
using DataFrames
using Statistics
using ProgressMeter
#using RCall

#=
ST-DBSCAN
@article{birant2007st,
  title={ST-DBSCAN: An algorithm for clustering spatial--temporal data},
  author={Birant, Derya and Kut, Alp},
  journal={Data \& knowledge engineering},
  volume={60},
  number={1},
  pages={208--221},
  year={2007},
  publisher={Elsevier}
}
https://www.sciencedirect.com/science/article/pii/S0169023X06000218

nb: output for `stnb()`
vals: Tuple (D, Δϵ)
  D: values Matrix
  Δϵ: difference between cluster mean and `D`
type: "GridCell", "Random"
  "GridCell": To narrow down the number of adjacent grid points, only points nb[n] with target (v[i] - v[n]) <= Δϵ are counted as adjacent points.
  "Random" : use nb[n] ## There are no regulations
minPts: minimum size of Neighbors. if Neighbors >= Δϵ then nb[i] will be a cluster.
=#
function st_dbscan(nb,
                   vals::Tuple...;
                   type::String = "GridCell",
                   minPts::Int64 = 10)
    
    ## set cluster number
    cluster::Int64 = 0;

    ## set cluster label array
    label = fill("", length(nb));

    ## ST-DBSCAN
    println("\nStart Clustering:  ", now())
    @showprogress for i in eachindex(nb)
        ## if label == "Noise" then not run.
        #println("Process: " * string(i) * "/" * string(length(nb)))
        if isempty(nb[i])
            #println("Skip: Number " * string(i) * " have not Neighbors.")
            label[i] = "Noise";
	elseif isempty(label[i]) || label[i] == "Noise"
            ## 格子データを用いる場合
            if type == "GridCell"
                # strClust = [];                
	        strClust =  map(vals) do x
	            nb[i][findall(<=(x[2]), abs.(x[1][i] .- x[1][nb[i]]))]
                end |> intersects;
                # strClust = searchNb(vals, nb[i], i);
            elseif type == "Random"
                strClust = nb[i];
            end
            ## minPtsよりクラスター数が少ない場合は "Noise" 判定
            if length(strClust) <= length(nb[i]) < minPts
                label[i] = "Noise";
            elseif minPts <= length(strClust) <= length(nb[i])
                ## up to cluster number
                cluster += 1;
                label[i] = string(cluster);
                println("\n --- Create Cluster : " * string(cluster));

                ## check cluster
                if type == "GridCell"
	            label[strClust] .= string(cluster);
                    label[setdiff(nb[i], strClust)] .= "Noise";
                elseif type == "Random"
                    label[nb[i]] .= string(cluster);
                end

                ## check node and Add Cluster.
                num = nb[i];
                while length(num) !== 0
                    nodeList = pmap(vals) do x

                        ## search cluster value
                        m = x[1][findall(==(string(cluster)), label)];

                        ## check nb =====
	                linkNb = unions(nb[num]);
                        if length(linkNb) == 1
	                    linkNb = linkNb[1]
                        end
                        if type == "GridCell"
                            fitNb = map(linkNb) do l
                                #nbs = searchNb(vals, nb[x], i);
                                #nbs = nb[x];
                                # nbs = Statistics.mean(vcat(m, nb[l]));
                                xd = x[1][nb[l]];
	                        μ = Statistics.mean(vcat(m, xd));
                                nbs = nb[l][findall(<=(x[2]), abs.(μ .- xd))];
                                if length(nbs) >= minPts
	                            return(nbs);
                                else
                                    return([]);
                                end
                            end |> unions;
                            if length(fitNb) == 1
	                        fitNb = fitNb[1]
                            end
                        elseif type == "Random"
                            fitNb = linkNb;
                        end
                        checkNb = setdiff(fitNb, num);
                        linkNb = nothing;
                        fitNb = nothing;

                        ## Fit Cluster:: search node =====
                        empNode = checkNb[isempty.(label[checkNb])];
                        # empNode = checkNb[findall(x -> isempty(x), label[checkNb]) ∩
                        #     findall(x -> x != "Noise", label[checkNb])];
                        
                        add = [];
                        for n in empNode
	                    μ = Statistics.mean(vcat(m, x[1][add]));
                            if abs(μ - x[1][n]) <= x[2]
	                        push!(add, n);
                            end
                            μ = nothing;
                        end
                        return(add)
                    end;
                    if length(nodeList) !== 0
	                node = intersects(nodeList);
                        nodeList = nothing;
                    else
                        node = [];
                    end
                    if length(node) !== 0
                        #println("- Add point for cluster-" * string(cluster) * " : " * join(node, " "))
                        label[node] .= string(cluster);
                        num = node;
                    else
                        num = [];
                    end
                end
            else
                label[i] = "Noise";
            end
        end
    end
    println("\n ", now(), " Completed.")
    return(label)
end
function st_dbscan(nb,
                   vals::Tuple;
                   type::String = "GridCell",
                   minPts::Int64 = 10)
    
    ## set cluster number
    cluster::Int64 = 0;

    ## set cluster label array
    label = fill("", length(nb));

    ## ST-DBSCAN
    println("\nStart Clustering:  ", now())
    @showprogress for i in eachindex(nb)
        ## if label == "Noise" then not run.
        #println("Process: " * string(i) * "/" * string(length(nb)))
        if isempty(nb[i])
            #println("Skip: Number " * string(i) * " have not Neighbors.")
            label[i] = "Noise";
	elseif isempty(label[i]) || label[i] == "Noise"
            ## 格子データを用いる場合
            if type == "GridCell"
                strClust = nb[i][findall(<=(vals[2]), abs.(vals[1][i] .- vals[1][nb[i]]))]
            elseif type == "Random"
                strClust = nb[i];
            end
            ## minPtsよりクラスター数が少ない場合は "Noise" 判定
            if length(strClust) <= length(nb[i]) < minPts
                label[i] = "Noise";
            elseif minPts <= length(strClust) <= length(nb[i])
                ## up to cluster number
                cluster += 1;
                label[i] = string(cluster);
                println("\n --- Create Cluster : " * string(cluster));

                ## check cluster
                if type == "GridCell"
	            label[strClust] .= string(cluster);
                    label[setdiff(nb[i], strClust)] .= "Noise";
                elseif type == "Random"
                    label[nb[i]] .= string(cluster);
                end

                ## check node and Add Cluster.
                num = nb[i];
                while length(num) !== 0
                    ## search cluster value
                    m = vals[1][findall(==(string(cluster)), label)];

                    ## check nb =====
	            linkNb = unions(nb[num]);
                    if length(linkNb) == 1
	                linkNb = linkNb[1]
                    end
                    if type == "GridCell"
                        fitNb = map(linkNb) do l
                            #nbs = searchNb(vals, nb[x], i);
                            #nbs = nb[x];
                            # nbs = Statistics.mean(vcat(m, nb[l]));
                            xd = vals[1][nb[l]];
	                    μ = Statistics.mean(vcat(m, xd));
                            nbs = nb[l][findall(<=(vals[2]), abs.(μ .- xd))];
                            if length(nbs) >= minPts
	                        return(nbs);
                            else
                                return([]);
                            end
                        end |> unions;
                        if length(fitNb) == 1
	                    fitNb = fitNb[1]
                        end
                    elseif type == "Random"
                        fitNb = linkNb;
                    end
                    checkNb = setdiff(fitNb, num);
                    linkNb = nothing;
                    fitNb = nothing;

                    ## Fit Cluster:: search node =====
                    empNode = checkNb[isempty.(label[checkNb])];

                    nodeList = [];
                    for n in empNode
	                μ = Statistics.mean(vcat(m, vals[1][nodeList]));
                        if abs(μ - vals[1][n]) <= vals[2]
	                    push!(nodeList, n);
                        end
                        μ = nothing;
                    end
                    if length(nodeList) !== 0
	                node = intersects(nodeList);
                        nodeList = nothing;
                    else
                        node = [];
                    end
                    if length(node) !== 0
                        #println("- Add point for cluster-" * string(cluster) * " : " * join(node, " "))
                        label[node] .= string(cluster);
                        num = node;
                    else
                        num = [];
                    end
                end
            else
                label[i] = "Noise";
            end
        end
    end
    println("\n ", now(), " Completed.")
    return(label)
end
# function st_dbscan(nb,
#                    vals::Tuple...;
#                    type::String = "GridCell",
#                    minPts::Int64 = 10)
    
#     ## set cluster number
#     cluster::Int64 = 0;

#     ## set cluster label array
#     label = fill("", length(nb));

#     ## ST-DBSCAN
#     println("\nStart Clustering:  ", now())
#     @showprogress for i in eachindex(nb)
#         ## if label == "Noise" then not run.
#         #println("Process: " * string(i) * "/" * string(length(nb)))
#         if isempty(nb[i])
#             #println("Skip: Number " * string(i) * " have not Neighbors.")
#             label[i] = "Noise";
# 	elseif isempty(label[i]) || label[i] == "Noise"
#             ## 格子データを用いる場合
#             if type == "GridCell"
#                 # strClust = [];                
# 	        strClust =  map(vals) do x
# 	            nb[i][findall(<=(x[2]), abs.(x[1][i] .- x[1][nb[i]]))]
#                 end |> intersects;
#                 # strClust = searchNb(vals, nb[i], i);
#             else
#                 strClust = nb[i];
#             end
#             ## minPtsよりクラスター数が少ない場合は "Noise" 判定
#             if length(strClust) <= length(nb[i]) < minPts
#                 label[i] = "Noise";
#             elseif minPts <= length(strClust) <= length(nb[i])
#                 ## up to cluster number
#                 cluster += 1;
#                 label[i] = string(cluster);
#                 println("\n --- Create Cluster : " * string(cluster));

#                 ## check cluster
#                 if type == "GridCell"
# 	            label[strClust] .= string(cluster);
#                     label[setdiff(nb[i], strClust)] .= "Noise";
#                 else
#                     label[nb[i]] .= string(cluster);
#                 end

#                 ## check node and Add Cluster.
#                 num = nb[i];
#                 while length(num) !== 0
#                     nodeList = map(vals) do x

#                         ## search cluster value
#                         m = x[1][findall(==(string(cluster)), label)];

#                         ## check nb =====
# 	                linkNb = unions(nb[num]);
#                         if length(linkNb) == 1
# 	                    linkNb = linkNb[1]
#                         end
#                         if type == "GridCell"
#                             fitNb = map(linkNb) do l
#                                 #nbs = searchNb(vals, nb[x], i);
#                                 #nbs = nb[x];
#                                 # nbs = Statistics.mean(vcat(m, nb[l]));
#                                 xd = x[1][nb[l]];
# 	                        μ = Statistics.mean(vcat(m, xd));
#                                 nbs = nb[l][findall(<=(x[2]), abs.(μ .- xd))];
#                                 if length(nbs) >= minPts
# 	                            return(nbs);
#                                 else
#                                     return([]);
#                                 end
#                             end |> unions;
#                             if length(fitNb) == 1
# 	                        fitNb = fitNb[1]
#                             end
#                         else
#                             fitNb = linkNb;
#                         end
#                         checkNb = setdiff(fitNb, num);
#                         linkNb = nothing;
#                         fitNb = nothing;

#                         ## Fit Cluster:: search node =====
#                         empNode = checkNb[isempty.(label[checkNb])];
#                         # if isempty(empNode)
# 	                #     print([]);
#                         # end

#                         add = [];
#                         for n in empNode
# 	                    μ = Statistics.mean(vcat(m, x[1][add]));
#                             if abs(μ - x[1][n]) <= x[2]
# 	                        push!(add, n);
#                             end
#                             μ = nothing;
#                         end
#                         return(add)
#                     end;
#                     if length(nodeList) !== 0
# 	                node = intersects(nodeList);
#                         nodeList = nothing;
#                     else
#                         node = [];
#                     end
#                     if length(node) !== 0
#                         #println("- Add point for cluster-" * string(cluster) * " : " * join(node, " "))
#                         label[node] .= string(cluster);
#                         num = node;
#                     else
#                         num = [];
#                     end
#                 end
#             else
#                 label[i] = "Noise";
#             end          
#         end
#     end
#     println("\n ", now(), " Completed.")
#     return(label)
# end
# @benchmark st_dbscan(nb, (v, Δϵ), minPts = 3)
# function fitClust(label, value, target, node, Δϵ)
#     node[findall(<=(Δϵ), abs.(target .- value[node]))]
# end



#=
地点別に隣接点行列を作成
x: x成分(経度)
y: y成分(緯度)
time: x, yに対応する時間情報ベクトル
eps0: 空間上の地点間距離の最小値
eps1: 空間上の地点間距離の最大値
eps2: 時間間隔の最大値
method: 距離の計算方法の指定
  "geo": 地理空間上の距離[km]として計算 use geoSailing
  "euclidian": ユークリッド距離として計算
=#
function stnb(x::Vector{Float64},
              y::Vector{Float64},
              time;
              eps0 = 0,
              eps1,
              eps2,
              method::String = "geo",
              type::String = "GridCell")
    
    # set output frame
    if type == "GridCell"
        len = length(x)
        timediff = map(time) do t
	    #abs.(time .- t);
            findall(x -> x <= eps2, abs.(time .- t))
        end
        ret = fill([], len, length(time))
        @sync Threads.@threads for i in eachindex(x)
	    g = findall(x -> eps0 <= x <= eps1, dist.(x, y, x[i], y[i], method = method))
            d = map(timediff) do t
                setdiff(vec(g .+ ((t .- 1) .* len)'), i)
            end
            ret[i, :] .= d
        end
    elseif type == "Random"
        ret = fill([], length(time))
        @sync Threads.@threads for i in eachindex(x)
	    ret[i] = findall(x -> eps0 <= x <= eps1, dist.(x, y, x[i], y[i], method = method)) ∩
                findall(x -> (time[i] - time[i]) < x <= eps2, abs.(time[i] .- time))
        end
    end

    return ret;
end
# function stnb(x::Vector{Float64},
#               y::Vector{Float64},
#               time;
#               eps0::Int64 = 0,
#               eps1::Int64,
#               eps2,
#               method::String = "geo")
    
#     # set output frame
#     len = length(x)
#     timediff = map(time) do t
# 	#abs.(time .- t);
#         findall(x -> x <= eps2, abs.(time .- t))
#     end
#     ret = fill([], len, length(time))
#     #@showprogress @inbounds @simd for i in eachindex(x)
#     @showprogress for i in eachindex(x)        
# 	g = findall(x -> eps0 <= x <= eps1, dist.(x, y, x[i], y[i], method = method))
#         d = map(timediff) do t
#             setdiff(vec(g .+ ((t .- 1) .* len)'), i)
#         end
#         ret[i, :] .= d
#     end
#     return ret;
# end
# @benchmark nb = stnb(x1, y1, time, eps1 = eps1, eps2 = eps2, method = "geo")


#=
地点間距離の算出
method:
  "euclidian": ユークリッド距離の計算
  "geo": 測地線航海算法距離の計算
=#
function dist(x1::Float64, y1::Float64,
              x2::Float64, y2::Float64;
              method::String="euclidian",
              x...)
    
    if method == "euclidian"
	√((x1-x2)^2 + (y1-y2)^2)
    elseif method == "geo"
        geoSailing(x1, y1, x2, y2; x...)
    end
end

# @benchmark dist.(x1, y1, x1', y1', method = "geo", Unit = "km", digit = 1)
# @benchmark for i in 1:length(x1)
#     dist.(x1, y1, x2[i], y2[i], method = "geo", Unit = "km", digit = 1)
# end

#=
測地線航海算法による回転楕円体上での直線距離計算アルゴリズム
https://www2.nc-toyama.ac.jp/WEB_Profile/mkawai/lecture/sailing/geodetic/geosail.html
=#
function geoSailing(Lon1::Float64, Lat1::Float64,
                    Lon2::Float64, Lat2::Float64;
                    A::Float64 = 6378137.0,
                    F::Float64 = 1/298.257222101,
                    Unit::String = "km",
                    digit::Int = 1)::Float64

    if Lon1 == Lon2 && Lat1 == Lat2
	return(0.0)
    end
    ## 扁平率
    #F = (A - B) / A
    B = A * (1 - F)

    ## ラジアン角度に変換
    lonA = deg2rad(Lon1)
    latA = deg2rad(Lat1)
    lonB = deg2rad(Lon2)
    latB = deg2rad(Lat2)

    ## 化成緯度に変換
    ϕA = atan(B / A * tan(latA))
    ϕB = atan(B / A * tan(latB))

    ## 球面上の距離
    χ = acos(sin(ϕA) * sin(ϕB) + cos(ϕA) * cos(ϕB) * cos(lonA - lonB))

    ## Lambert-Andoyer補正
    Δρ = F / 8 * ((sin(χ) - χ) * (sin(ϕA) + sin(ϕB)) ^ 2 / cos(χ / 2) ^ 2 - (sin(χ) + χ) * (sin(ϕA) - sin(ϕB)) ^ 2 / sin(χ / 2) ^ 2)
    if isnan(Δρ)
	Δρ = 0.0
    end
    
    ## 測地線長
    ρ = A * (χ + Δρ)

    if Unit == "km"
        return(round(ρ / 1000, digits = digit))
    else
        return(round(ρ, digits = digit))
    end        
end
# geoSailing.(x1, y1, x2, y2, Unit = "km", digit = 0)
# geoSailing.(x2, y2, x2', y2', Unit = "m", digit = 1)

            
#=
条件に合致する隣接点をピックアップ
=#
function searchNb(vals, nb, i)
    clList =  map(vals) do x
	nb[findall(<=(x[2]), abs.(x[1][i] .- x[1][nb]))]
    end
    intersects(clList);
end

function merge2(x::Tuple)
    if length(x) == 1
        x
    elseif length(x) == 2
        vcat(x[1], x[2])
    else
        vcat(x[begin], merge2(x[begin+1:end]))
    end
end
function intersects(x)
    if length(x) <= 1
	x
    elseif length(x) == 2
        intersect(x[1], x[2])
    else
        intersect(x[1], intersects(x[2:end]))
    end
end
function unions(x)
    if length(x) <= 1
	x
    elseif length(x) == 2
        union(x[1], x[2])
    else
        union(x[1], unions(x[2:end]))
    end
end
function rbind(x)
    if length(x) == 2
	vcat(x[1], x[2])
    else
        vcat(x[1], rbind(x[2:end]))
    end
end

function getSkipData(x, y;
                     skip_x::Int64 = 0,
                     skip_y::Int64 = 0,
                     interval_x::Int64 = 5,
                     interval_y::Int64 = 5)
    xskip = sort(unique(x), rev = false);
    xs = xskip[1 + skip_x:interval_x:length(xskip)];
    yskip = sort(unique(y), rev = true);
    ys = yskip[1 + skip_y:interval_y:length(yskip)];
    boolxy = findall(x -> x ∈ xs, x) ∩ findall(y -> y ∈ ys, y)
    # zz = map(z) do zz
    #     zz[boolxy]
    # end
    geo = allcombinations(DataFrame,
                          x = xs,
                          y = ys)

    #return (x = geo.x, y = geo.y, z = zz)
    return (df = geo, num = boolxy)
end
