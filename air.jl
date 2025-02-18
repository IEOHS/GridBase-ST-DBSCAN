using Dates: format
using RCall: CategoricalArrays
using DataFrames: table_transformation
using ArchGDAL: assingle
using Base: stop_reading, fill_to_length
using BenchmarkTools
using Dates
using DataFrames
using GeoDataFrames
using Shapefile
using Statistics
using NetCDF
using Distributed
using GeoFormatTypes
using Statistics
using GDAL
using ArchGDAL


"""
風速の計算
"""
function ws(u::Float64, v::Float64; digits::Int64 = 2)
    return round(sqrt(u ^ 2 + v ^ 2); digits = digits);
end

"""
風向の計算
"""
function wd(u::Float64, v::Float64; digits::Int64 = 2)
    deg = atan(v / u) * 180 / π;
    return round(ifelse(u >= 0, 270 - deg, 90 + deg); digits = digits);
end


GC.gc()


#=
スキップデータを作成
=#
function getSkipedAddress(geo::Matrix;
                          skip_x::Int64 = 4,
                          skip_y::Int64 = 4,
                          bbox::NamedTuple = (xmin::Float64, ymin::Float64,
                                              xmax::Float64, ymax::Float64));
  
    ## get intersects
    if 0 <= length(bbox) < 4 && 4 < length(bbox)
        println("skip clip process.")
        clip_area = [CartesianIndex((x, y)) for x in 1:size(geo, 1), y in 1:size(geo, 2)];
    elseif length(bbox) == 4
        boundarybox = GeoDataFrames.createpolygon([(bbox.xmin, bbox.ymin)
                                                   (bbox.xmin, bbox.ymax)
                                                   (bbox.xmax, bbox.ymax)
                                                   (bbox.xmax, bbox.ymin)
                                                   (bbox.xmin, bbox.ymin)]);
        clip_area = findall(x -> !ArchGDAL.isempty(x),
                            GeoDataFrames.intersection.([boundarybox], geo));        
    end
    
    ## get skip data
    if skip_x != 1 && skip_y != 1
        skips = [CartesianIndex((x, y)) for x in 1:skip_x:size(geo, 1), y in 1:skip_y:size(geo, 2)];
        drop_address = [[CartesianIndex((xx, yy)) for xx in x:(x + skip_x - 1), yy in y:(y + skip_y - 1)]
                        for x in 1:skip_x:size(geo, 1), y in 1:skip_y:size(geo, 2)];
        
        num = findall(x -> x in clip_area, skips);
        return (address = skips[num],
                same_address = drop_address[num]);
    else
        return (address = clip_area,
                same_address = nothing);
    end
end


#=
Netcdfファイルのgpvデータを変換する
=#
function gpv2geo(geo)
    poly = [];
    for i in 1:(size(geo, 1) - 1), j in 1:(size(geo, 2) - 1)
        push!(poly, [geo[i:(i + 1), j:(j + 1)][[1, 2, 4, 3, 1]]])
    end
    return reshape(createpolygon.(poly), (size(geo, 2) - 1, size(geo, 1) - 1)) |> permutedims
end 
function geocentroid(geo)
    poly = [];
    for i in 1:(size(geo, 1) - 1), j in 1:(size(geo, 2) - 1)
        push!(poly, (mean(map(x -> x[1], geo[i:(i + 1), j:(j + 1)])),
                     mean(map(x -> x[2], geo[i:(i + 1), j:(j + 1)]))))
    end
    return reshape(poly, (size(geo, 2) - 1, size(geo, 1) - 1)) |> permutedims
end
@generated function gpv2var(vars; fun = mean)
    ##:(fun(vars))
    quote
        meanarray = [];   
	for i in 1:(size(vars, 1) - 1), j in 1:(size(vars, 2) - 1)
            push!(meanarray, fun(vars[i:(i + 1), j:(j + 1)]))
        end
        reshape(meanarray, (size(vars, 2) - 1, size(vars, 1) - 1)) |> permutedims
    end
end


#=
正規表現でファイルを絞り込む関数
=#
function searchfile(path::AbstractString, keywords::Regex)
    filter(x -> occursin(keywords, x),
           readdir(exhomedir(path), join = true))
end
function exhomedir(path::AbstractString)
    if occursin("~", path)
	return replace(path, "~" => homedir())
    else
        return path
    end
end

#=
時刻を取得
=#
function ncgettime(ncfile::AbstractString)
    nc_hour = ncvarget(ncfile, "time");
    st_time = DateTime(SubString(replace(NetCDF.ncgetatt(ncfile,
                                                         "time",
                                                         "units"),
                                         "hours since " => ""),
                                 1, 10));
    return st_time .+ Hour.(nc_hour)
end

#=
netCDFファイルから情報を取得する関数
=#
function ncvarget(ncfile::AbstractString,
                  var::AbstractString)
    val = NetCDF.ncread(ncfile, var);
    attr = map(["scale_factor", "add_offset"]) do x
        try
	    NetCDF.ncgetatt(ncfile, var, x);
        catch
        end
    end;
    if !any(@. isnothing(attr))
	return (val .* attr[1]) .+ attr[2];
    else
        return val;
    end
end



"""
Check airmass
lower: lower size for airmass[km]
upper: upper size for airmass[km]
"""
function airmass(lon, lat, clusterLabel;
                 lower::Int64 = Int64(round(√(100 ^ 2 + 100 ^ 2), digits = 0)),
                 upper::Int64 = 4000)
    cl = setdiff(unique(clusterLabel), ["Noise"]);
    println("----------------------------------------")
    println("check AirMass.  $(now())")
    ret = fill("Noise", length(clusterLabel));
    cluster::Int64 = 0;
    for lab in cl
	n = findall(x -> x == lab, clusterLabel);
        bbox = [minimum(lon[n]),
                minimum(lat[n]),
                maximum(lon[n]),
                maximum(lat[n])];
        clsize = geoSailing(bbox[1], bbox[2], bbox[3], bbox[4]);
        if lower <= clsize <= upper
            cluster += 1;
            ret[n] .= "airmass_$(string(cluster))";
	    println("$(round(lower, digits = 1)) <= Cluster: $lab <= $upper [km]. Create Airmass = airmass_$cluster")
        else
            println("Cluster: $lab is under $lower [km] or grater than $upper [km].  Delete cluster.")
        end
    end
    return ret
end
    
"""
温位θの計算
t: 気温 °C
p: 気圧 hPa
p₀: 基本気圧 hPa
Rd:: 乾燥空気の気体定数 J/K⋅kg
Cₚ: 定圧比熱 J/K⋅kg

Unit: K
"""
function PotentialTemp(t::Float64 = 20.0,
                       p::Float64 = 1013.25,
                       p₀::Float64 = 1000.0,
                       Rd:: Int64 = 287,
                       Cₚ::Int64 = 1004)
    (t + 273.15) * (p / p₀) ^ (-Rd / Cₚ)
end


"""
比湿qの計算
mᵣ: 混合比 g/kg

Unit: g/kg
"""
function SpecificHum(mᵣ::Float64)
     mᵣ / (1 + mᵣ);
end


"""
混合比mᵣの計算
t: 気温 °C
RH: 相対湿度 %
p: 気圧 hPa
eₛ: 飽和水蒸気圧 hPa

Unit: g/kg
"""
function MixRatio(RH::Float64 = 50.0,
                  p::Float64 = 1013.25,
                  eₛ::Float64 = 0.0)
    e = eₛ * RH / 100.0;
    622 * e / (p - e) ## g/kg
end

"""
飽和水蒸気圧eₛの計算
t: 気温 °C
Tₖ: 沸点 K
eₛₜ: 基準気圧 hPa

Unit: hPa
"""
function GoffGratch(t::Float64 = 20.0,
                    Tₖ::Float64 = 373.16,
                    eₛₜ::Float64 = 1013.25)
    T = t + 273.15;
    e = -7.90298(Tₖ / T - 1) +
        5.02808log10(Tₖ / T) -
        1.3816 * 10 ^ (-7) * (10 ^ (11.344 * (1 - T / Tₖ)) - 1) +
        8.1328 * 10 ^ (-3) * (10 ^ (-3.49149 * (Tₖ / T - 1)) - 1) +
        log10(eₛₜ)
    return(10 ^ e)
end


"""
飽和水蒸気量の計算
t: 気温 °C
RH: 相対湿度 %
eₛ: 飽和水蒸気圧 hPa

Unit: g/m³
"""
function WaterVapor(t::Float64 = 20.0,
                    RH::Float64 = 50.0,
                    eₛ::Float64 = 50.0)
    ## 飽和水蒸気量
    (217.0 * eₛ / (t + 273.15)) * RH / 100.0
end
