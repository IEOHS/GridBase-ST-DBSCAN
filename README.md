ST-DBSCAN for Julia and R
=========================

ST-DBSCAN[1]は密度ベースのクラスタリング手法であり、時空間データの解像度はそのままにクラスタリングを行い、ホットスポット検出など時空間データの解析に広く利用されている手法です。  
このプロジェクトでは、Julia言語[2]及びR言語[3]でST-DBSCANを実装し、データ解析に利用しやすくすることを目的としています。  

ST-DBSCANのアルゴリズム
=======================

ST-SBSCAN法のアルゴリズム[4]は引用文献を元に記述しました。また、このプロジェクトではグリッドデータにも対応できるよう、アルゴリズムを一部修正しています。  
具体的には、クラスター認定されるためには地点 *X* の隣接点 *Y* の数が
MinPts以上であることが条件 ( \|*Y*\| ≧ *M**i**n**P**t**s* )
となっていますが、グリッドデータは常に一定数の隣接点を持つことからこの条件では最適な隣接点を求めることができないと思われます。  
そこで、 地点 *X* と隣接点 *Y* の値 *D* の差が *Δ**ϵ* 以上となる条件
*Y*<sub>*ϵ*</sub> = (*X* − *Y*) ≧ *Δ**ϵ* の数が MinPts 以上となる条件 (
\|*Y*<sub>*ϵ*</sub>\| ≧ *m**i**n**P**t**s*
)に変更してクラスター判定を行うよう修正しました。

Demo
====

関数のインポート
----------------

関数定義ファイルをインポートしてください。

``` julia
include("./st_dbscan.jl")
```

``` r
source("./st_dbscan.R")
```

ランダム地点データ
------------------

### 計算用データの作成

1.  for Julia

    ``` julia
    using Dates
    using DataFrames
    x = Float64.(round.(rand(100) * 100, digits = 1));
    y = Float64.(round.(rand(100) * 100, digits = 1));
    t = DateTime("2024-01-01T00:00:00") .+ Hour.(0:6:23);
    geo = DataFrame(x = x, y = y, t = repeat(t, 25));
    geo.D = rand(nrow(geo));
    first(geo, 10)
    ```

2.  for R

    ``` r
    x <- round(abs(runif(100)) * 100, digits = 1)
    y <- round(abs(runif(100)) * 100, digits = 1)
    t <- as.POSIXct("2024-01-01 00:00:00", tz = "JST") + 3600 * seq(0, 23, by = 6)
    geo <- data.frame(x = x, y = y, t = rep(t, 25))
    geo$D <- abs(runif(100))
    head(geo) 
    ```

### 隣接点情報の作成

隣接点情報は `stnb()` 関数を使って計算します。
グリッドデータとランダム地点データでは `stnb()`
の引数に指定するデータ数が異なりますので注意してください。  

1.  for Julia

    ``` julia
    nb = stnb(geo.x, geo.y, geo.t,
              eps1 = 20, eps2 = Hour(6), method = "euclidian", type = "Random")
    ```

2.  for R

    ``` r
    nb <- stnb(geo$x, geo$y, geo$t,
               eps1 = 20, eps2 = 3600*6, method = "euclidian", type = "Random")
    ```

### ST-DBSCANの実行

クラスタリングの結果は文字列ベクトルで出力されます。

1.  for Julia

    ``` julia
    geo.clust = st_dbscan(nb, (geo.D, 0.3); minPts = 5, type = "Random")
    first(geo, 10)
    ```

    ``` julia
    10×5 DataFrame
     Row │ x        y        t                    D           clust  
         │ Float64  Float64  DateTime             Float64     String 
       1 │     2.6     14.7  2024-01-01T00:00:00  0.0954608   Noise──
       2 │    50.3     12.2  2024-01-01T06:00:00  0.303733    1
       3 │    65.5     11.5  2024-01-01T12:00:00  0.102751    10
       4 │    20.1     50.4  2024-01-01T18:00:00  0.675904    2
       5 │    14.2      3.3  2024-01-01T00:00:00  0.702335    Noise
       6 │    17.8     66.7  2024-01-01T06:00:00  0.299661    3
       7 │    21.4     31.4  2024-01-01T12:00:00  0.192165    2
       8 │    71.2     66.5  2024-01-01T18:00:00  0.111099    Noise
       9 │    44.9     12.5  2024-01-01T00:00:00  0.8918      1
      10 │    66.4     56.1  2024-01-01T06:00:00  0.00891727  4
    ```

2.  for R

    ``` r
    geo$clust <- st_dbscan(nb = nb,
                           vals = list(list(D = geo$D, delta_eps = 0.1)),
                           type = "Random", minPts = 5)
    ```

グリッドデータ
--------------

### 計算用データの作成

1.  for Julia

    ``` julia
    using DataFrames
    x = Float64.(130.0:1:140.0);
    y = Float64.(30.0:1:40);
    t = collect(1:4)
    geo = allcombinations(DataFrame, x = x, y = y)
    first(geo, 10)
    ```

2.  for R

    ``` r
    x <- seq(130, 140, by = 1)
    y <- seq(30, 40, by = 1)
    t <- 1:4
    geo <- with(expand.grid(x, y), {
      data.frame(x = Var1, y = Var2)
    })
    ```

### 隣接点情報の作成

1.  for Julia

    ``` julia
    nb = stnb(geo.x, geo.y, t,
              eps1 = 144, eps2 = 1, method = "geo", type = "GridCell")
    ## 時間まで含めてデータブレームを再生成
    geo = allcombinations(DataFrame, x = x, y = y, t = t)
    geo.D = rand(nrow(geo)) * 100;
    first(geo, 10)
    ```

    ``` julia
    10×5 DataFrame
     Row │ x        y        t      D        
         │ Float64  Float64  Int64  Float64  
    ─────┼───────────────
       1 │   130.0     30.0      1  19.1332  
       2 │   131.0     30.0      1  80.7502  
       3 │   132.0     30.0      1  33.5397  
       4 │   133.0     30.0      1  59.749   
       5 │   134.0     30.0      1   1.76694 
       6 │   135.0     30.0      1  89.2776  
       7 │   136.0     30.0      1  53.4877  
       8 │   137.0     30.0      1  14.6705  
       9 │   138.0     30.0      1  11.5225  
      10 │   139.0     30.0      1  74.7262  
    ```

2.  for R

    ``` r
    nb <- stnb(geo$x, geo$y, t,
               eps1 = 144, eps2 = 1, method = "geo", type = "GridCell")
    geo <- with(expand.grid(x, y, t), {
      data.frame(x = Var1, y = Var2, t = Var3)
    })
    geo$D <- abs(runif(nrow(geo))) * 100
    ```

### ST-DBSCANの実行

1.  for Julia

    ``` julia
    geo.clust = st_dbscan(nb, (geo.D, 20); minPts = 11, type = "GridCell")
    first(geo, 10)
    ```

    ``` julia
    10×5 DataFrame
     Row │ x        y        t      D         clust  
         │ Float64  Float64  Int64  Float64   String 
    ─────┼───────────────────────────────────────────
       1 │   130.0     30.0      1  19.1332   Noise
       2 │   131.0     30.0      1  80.7502   Noise
       3 │   132.0     30.0      1  33.5397   Noise
       4 │   133.0     30.0      1  59.749    Noise
       5 │   134.0     30.0      1   1.76694  Noise
       6 │   135.0     30.0      1  89.2776   Noise
       7 │   136.0     30.0      1  53.4877   Noise
       8 │   137.0     30.0      1  14.6705   Noise
       9 │   138.0     30.0      1  11.5225   Noise
      10 │   139.0     30.0      1  74.7262   Noise
    ```

2.  for R

    ``` r
    geo$clust <- st_dbscan(nb = nb,
                           vals = list(list(D = geo$D, delta_eps = 20)),
                           type = "GridCell", minPts = 11)
    ```

Licence
=======

このスクリプトファイルは MIT ライセンスで公開します。

[1] BIRANT, Derya; KUT, Alp. ST-DBSCAN: An algorithm for clustering
spatial–temporal data. Data & knowledge engineering, 2007, 60.1:
208-221.
<https://www.sciencedirect.com/science/article/pii/S0169023X06000218>

[2] <https://julialang.org/>

[3] <https://cran.r-project.org/>

[4] BIRANT, Derya; KUT, Alp. ST-DBSCAN: An algorithm for clustering
spatial–temporal data. Data & knowledge engineering, 2007, 60.1:
208-221.
<https://www.sciencedirect.com/science/article/pii/S0169023X06000218>
