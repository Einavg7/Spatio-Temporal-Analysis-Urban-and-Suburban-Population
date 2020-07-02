---
title: "Urban and Suburban Population Density Changes between 2010-2018 in Massachusetts, USA"
author: "Einav Grinberg"
date: "February 26, 2020"
output: 
  html_document: 
    keep_md: true
---



***

## Introduction

The United States of America has gone through a process of change regarding human settlement development throughout the years. 
The mobility of population in the USA can be divided into three periods. In the begining, people arrived to the eastern coastal area, and started migrating to the west, creating a dispersing population all over the country. Subsequently, during the 19th century and mainly in the 1950's population distribution shifted, and more people moved to rural and suburban areas and avoided living in big cities. The third period starts with the rise of the New Urbanism movement that was formed in the USA in the early 1980's and continues to influence human settlement patterns in present times (Federal Research Division of the Library of Congress, n.d.). 
In a Report published by the United Nations Population Division in 2019, it is stated that more than half of the world's people are living in urban areas. The  world's  urban population  is  now  close  to 4.2 billion and is expected to reach 6.7 billion in 2050. When examining urban population trends in the USA, in 1790 less than 10% of the population is defined as urban compared to 2018, where 80% of the population is estimated as urban. By 2050, the projected percentage of the urban population in the USA is estimated to be 90%. Urban population growth is divided in to two forms of urbanization. The first form is described as high dense urban area, meaning less urban sprawl. The second form is characterized by low dense urban area, meaning more urban sprawl. These changes have a great impact on land-use policies. The efforts to conserve open space are attempting to address the competition between urban, agriculture and nature conservation while population is continuing to increase and the land area remains the same (Paige, Ryan, Lerman & Tooke, 2011). 
Massachusetts is composed by big main cities such as Boston, Springfield, Worcester, Lowell and Cambridge and many suburban and rural settlements. As a result of Massachusetts's geographic location on the east coast of the USA many immigrants have entered and settled in different areas throughout the years. In addition to the aforementioned changes in population distribution in the USA, the most populated cities in Massachusetts have a reputation for a highly dense urban population core, while, adjacent towns and suburbs are much less dense, and urban sprawl is becoming an urging issue for land-use policy makers in the state of Massachusetts (Cox, 2015). 
This study's objective is to further analyze population density changes in Massachusttets between the years 2010-2018. The first hypothesis of this study is that cities population density is increasing compared to population density in towns and town-cities between 2010-2018. The second hypothesis of this study is that  population density in main cities is increasing compared to population density in their neighbour towns and town-cities between 2010-2018.


## Methodology 

Firstly, population change in all of Massachusetts was temporaly analyzed between 2010-2018. Next the population was divided by settlement type; city, town and town-city. 
The types are defined by Massachusetts Secretary of State Office, town-city is a town with city form of government.
The population growth rate between these years was calculated using $\frac{pop_t-pop_{t-1}}{pop_{t-1}}*{100}$ for each settlement type. After adding the spatial data, population density was calculated for each settlement using $\frac{pop_t}{area}$ in $\frac{1}{km^2}$.
The density changes have been visualized and temporaly analyzed between 2010-2018 for all of Massachusetts and for settlement types. To test for spatial autocorrelation between all the settlements in Massachusetts two types of neighbours and spatial weights have been used. The first set of neighbours was created using queen contiguous neighbours. This method is calculated by each polygon centroid checking its distance to at least one other polygon centroid. 
The second set of neighbours were calculated using ${k}$-nearest neighbours using the equation ${k}=\sqrt{n}$ where ${n}$ is the number of observations, as a threshold value. Next, the spatial weights were calculated using binary style that gives a weight of unity to each neighbour relationship, and upweights units with no boundaries on the edge of the study area. The spatial weights were calculated without Massachusetts islands using the zero.policy argument. After creating two sets of spatial weight based on the different neighbour calculations, Moran's I test was calculated for every year between all Massachusetts settlements. For further analysis, the five most populated cities and there neighbours were selected. First, a temporal analysis of density population changes was examined. As mentioned in the paper Spatial Regression Models for Demographic Analysis by Chi & Zhu (2007) the distance and contiguous based spatial weight matrix tend to be less compitable regarding census units, because they usually make too many neighbors in urban areas and too few neighbors in rural areas. Therefore, ${k}$-nearest neighbours based spatial weights were used for the main cities and their neighbours, the spatial weights were calculated using the equation ${k}=\sqrt{n}$ as a threshold value. With the spatial weights created, Moran's I test was computed to analyze spatial autocorrelation between the main cities and their neighbours population density for every year from 2010 to 2018. 

## Data

* The Population data was downloaded as a csv file from UMass State Data Center: http://www.donahue.umassp.edu/business-groups/economic-public-policy-research/massachusetts-population-estimates-program.

* The spatial data was downloadad as a shapefile from the Massachusetts Government Data Portal: https://docs.digital.mass.gov/dataset/massgis-data-county-boundaries.





### Population Data

The population data contains a table of 351 cities, towns and town-cities and their population numbers between the years 2010-2018.

When plotting the time series of the population numbers, the plot displays an increase in the total number of population in Massachusetts between 2010-2018.

![](FinalAssignment_files/figure-html/total population-1.png)<!-- -->



```r
pop_data3 = pop_data2
pop_data3[,352] = NULL
pop_data3 = rbind(type, pop_data3)
towns = pop_data3[, grepl("^t", pop_data3[1,])]
towns = towns[2:10,]
towns = as.data.frame(lapply(towns, as.numeric))
towns$totalPopulation = rowSums(towns)
towns = cbind(years, towns)
town_cities = pop_data3[, grepl("^T", pop_data3[1,])]
town_cities = town_cities[2:10,]
town_cities = as.data.frame(lapply(town_cities, as.numeric))
town_cities$totalPopulation = rowSums(town_cities)
town_cities = cbind(years, town_cities)
cities = pop_data3[, grepl("^c", pop_data3[1,])]
cities = cities[2:10,]
cities = as.data.frame(lapply(cities, as.numeric))
cities$totalPopulation = rowSums(cities)
cities = cbind(years,cities)
towns_ts =  ts(towns, start = 2010, end = 2018, frequency = 1)
town_cities_ts = ts(town_cities, start = 2010, end = 2018, frequency = 1)
cities_ts = ts(cities, start = 2010, end = 2018, frequency = 1)

#calculate total growth rate for each type
library(tis)
```

#### Annual Population Growth Rate in Massachusetts Towns

The annual population growth rate in Massachusetts towns between 2010 to 2011 is ~0.7% and then declines by ~0.07% and by 2012, increases to the highest percent of population growth observed; ~0.73%. From 2012 to 2015 there is a large decrease in population growth. Between 2015 and 2016 there is an increase again that remains stable at ~0.62% until 2018.


```r
plot(growth.rate(towns_ts[,299]), type = 'l', main = paste("Annual Population Growth Rate in Massachusetts Towns", "\nbetween 2010-2018"), ylab = "Growth Rate Percentage")
```

![](FinalAssignment_files/figure-html/growth-1.png)<!-- -->

#### Annual Population Growth Rate in Massachusetts Town-Cities

The annual population growth rate in Massachusetts town-cities between 2010 to 2014 declines by ~0.07% and then increases slightly. From 2016 to 2018 there is a large decrease in population growth that eventually decrease to ~0.02%.


```r
plot(growth.rate(town_cities_ts[,15]), type = 'l', main = paste("Annual Population Growth Rate in Massachusetts Town-Cities", "\nbetween 2010-2018"), ylab = "Growth Rate Percentage")
```

![](FinalAssignment_files/figure-html/growth town city-1.png)<!-- -->

#### Annual Population Growth Rate in Massachusetts Cities

The annual population growth rate in Massachusetts cities between 2010 to 2011 is ~0.7%. Between 2011 to 2014 the population increases by ~0.07%. From 2014 to 2015 there is a large decrease in population growth. Between 2015 to 2017 there is another slight decrease by ~0.05%, that eventually increases back to 0.5% up til 2018. 


```r
plot(growth.rate(cities_ts[,42]), type = 'l', main = paste("Annual Population Growth Rate in Massachusetts Cities", "\nbetween 2010-2018"), ylab = "Growth Rate Percentage")
```

![](FinalAssignment_files/figure-html/growth city, warning-1.png)<!-- -->

### Spatial Data

The spatial data is a shapefile with the polygons of 351 cities, towns and town-cities that was created using the Municipal Boundaries from Census 2010.

The first plot shows the distribution of Massachusetts towns, cities and town-cities. The second plot displays the population change in each polygon between 2010-2018 using the Jenks natural breaks classification method, which reduces the variance within classes and maximizes the variance between classes. 


```r
#read polygon shapefile
mass = read_sf("data/CENSUS2010TOWNS_SHP/FIXED_MASS_POLY.shp", stringsAsFactors = F)
mass = st_transform(mass, crs = 4326)

#prepare for join pop data and polygons
pop_data = as.data.frame(lapply(pop_data, gsub, pattern = ",", replacement = ""), stringsAsFactors = F)
pop_data[c(1,4:12)] <- sapply(pop_data[c(1,4:12)],as.numeric)
names(pop_data)[1] = "TOWN_ID"
mass_pop = left_join(mass, pop_data, by = "TOWN_ID")

#plot poly by type
plot(mass_pop["type"], key.pos = 1, main = "Massachusetts Towns, Town-Cities and Cities")
```

![](FinalAssignment_files/figure-html/polygons-1.png)<!-- -->

```r
#plot population 2010-2018
plot(mass_pop[24:32], key.pos = 1, breaks = "jenks")
```

![](FinalAssignment_files/figure-html/polygons-2.png)<!-- -->


## Population Density Analysis

Population Density plotted for each polygon for every year and classified using jenks breaks.


```r
#calculate population density for every year
library(units)
mass_pop$Area_Km = st_area(mass_pop)
mass_pop$Area_Km = set_units(mass_pop$Area_Km, "km^2")
mass_pop2 = mass_pop %>% mutate(density_2010 = X2010/Area_Km, density_2011 = X2011/Area_Km, density_2012 = X2012/Area_Km, density_2013 = X2013/Area_Km, density_2014 = X2014/Area_Km, density_2015 = X2015/Area_Km, density_2016 = X2016/Area_Km, density_2017 = X2017/Area_Km, density_2018 = X2018/Area_Km)

plot(mass_pop2[34:42], key.pos = 1, breaks = "fisher")
```

![](FinalAssignment_files/figure-html/density data-1.png)<!-- -->

```r
mass_density = mass_pop2[,c(22:23, 34:42)]
mass_density2 = as.data.frame(t(mass_density))
names(mass_density2) = mynames
mass_density2 = mass_density2[2:11,]
mass_density2 = mass_density2[2:10,]
mass_density2 = as.data.frame(lapply(mass_density2, as.numeric))
mass_density2$total_density = rowSums(mass_density2)
mass_density2 = cbind(years, mass_density2)
```

The population density between 2010 and 2018 in Massachusetts is increasing. In 2010 the population density per square kilometer was 171618.3 and in 2018 the population density per square kilometer was 181478.5.

The autocorrelation plot for Massachusetts population density displays that in lag 1 (year - 1) there is a correlation between density population of ~0.6 regarding the observed year. After this year the correlation in lag 2 is not significant. 


```r
#population density trend
plot(years, mass_density2[,353], ylab = "Population Density", main = "Population Density in Massachusetts between 2010-2018", xlab = "Years", type = 'l')
```

![](FinalAssignment_files/figure-html/pop_dense trend and acf-1.png)<!-- -->

```r
#acf for density
acf(mass_density2[,353], main = "Autocorrelation Massachusetts Population Density")
```

![](FinalAssignment_files/figure-html/pop_dense trend and acf-2.png)<!-- -->

## Population Density Analysis by Settlement Type


```r
mass_density3 =  as.data.frame(t(mass_density))
mynamesdensity = mass_density3[1,]
names(mass_density3) = unlist(mynamesdensity)
mass_density3 = mass_density3[2:11,]

towns_density = mass_density3[, grepl("town", mass_density3[1,])]
towns_density = towns_density[2:10,]
towns_density = as.data.frame(lapply(towns_density, as.numeric))
towns_density$total_density = rowSums(towns_density)
towns_density = cbind(years,towns_density)

town_cities_density = mass_density3[, grepl("Town city", mass_density3[1,])]
town_cities_density = town_cities_density[2:10,]
town_cities_density = as.data.frame(lapply(town_cities_density, as.numeric))
town_cities_density$total_density = rowSums(town_cities_density)
town_cities_density = cbind(years,town_cities_density)

cities_density = mass_density3[, grepl("\\Scity", mass_density3[1,])]
cities_density = cities_density[2:10,]
cities_density = as.data.frame(lapply(cities_density, as.numeric))
cities_density$total_density = rowSums(cities_density)
cities_density = cbind(years,cities_density)


townsdens_ts =  ts(towns_density, start = 2010, end = 2018, frequency = 1)
town_citiesdens_ts = ts(town_cities_density, start = 2010, end = 2018, frequency = 1)
citiesdens_ts = ts(cities_density, start = 2010, end = 2018, frequency = 1)
```

#### Annual Population Density Growth Rate in Massachusetts Towns

The annual population density growth rate in Massachusetts towns between 2010 to 2011 is ~0.72%. Between 2011 to 2012 the density declines slightly, and then increases to the highest percent of density growth observed; ~0.78%. From 2014 to 2016 there is a large decrease in density growth. Between 2016 and 2018 there is an increase again up til 2018 that is ~0.58%.


```r
plot(growth.rate(townsdens_ts[,299]), type = 'l', main = paste("Annual Population Density Growth Rate in Massachusetts Towns", "\n between 2010-2018"), ylab = "Growth Rate Percentage")
```

![](FinalAssignment_files/figure-html/townsdens-1.png)<!-- -->

#### Annual Population Density Growth Rate in Massachusetts Town-Cities

The annual population density growth rate in Massachusetts town-cities between 2010 to 2011 is ~0.86%. Between 2011 to 2012 the density increase by ~0.03%, and then by 2013 declines to 0.6% of density growth. From 2013 to 2014 there is a large increase in density growth to ~1.3%. Between 2014 and 2015 there is a significant decrease by ~1%. Between 2015 and 2016 the density is increasing up til 2017 when the density decreases by ~0.03% until 2018.


```r
plot(growth.rate(town_citiesdens_ts[,16]), type = 'l', main = paste("Annual Population Density Growth Rate in Massachusetts Town-Cities", "\n between 2010-2018"), ylab = "Growth Rate Percentage")
```

![](FinalAssignment_files/figure-html/townscitiesdens-1.png)<!-- -->

#### Annual Population Density Growth Rate in Massachusetts Cities

The annual population density growth rate in Massachusetts cities between 2010 to 2011 is ~0.89%. Between 2011 to 2013 the density increase by ~0.18%, and then by 2017 the density decreases to ~0.37%. From 2017 to 2018 there is a slight increase in density growth to ~0.41%.


```r
plot(growth.rate(citiesdens_ts[,42]), type = 'l', main = paste("Annual Population Density Growth Rate in Massachusetts Cities", "\n between 2010-2018"), ylab = "Growth Rate Percentage")
```

![](FinalAssignment_files/figure-html/citiesdens-1.png)<!-- -->


## Spatial Autocorrelation for Massachusetts Population Density 2010-2018

#### Spatial autocorrelation for contiguous neighbours - Massachusetts Density

The contiguties plot describes the calculated contiguous queen neighbours from each Massachusetts polygon centroid.

The contiguous neighbours are calculated using the queen method. Next, the spatial weights are calculated using zero.policy to disclude Massachusetts islands.  


```r
library(spdep)

#Contiguous queen neighbours
nb_q <- poly2nb(mass_density, queen=TRUE)
nb_q
```

```
## Neighbour list object:
## Number of regions: 351 
## Number of nonzero links: 1810 
## Percentage nonzero weights: 1.469144 
## Average number of links: 5.156695 
## 2 regions with no links:
## 118 179
```

```r
plot(st_geometry(mass_density), border = 'grey', main = "Massachusetts Queen-style Population Density Contiguities")
plot(nb_q, st_centroid(st_geometry(mass_density)), pch = 3, cex = .2, add = TRUE, randomisation = FALSE)
```

![](FinalAssignment_files/figure-html/spatial autocorrelation contiguous-1.png)<!-- -->

```r
#spatial weights for contiguous neighbours
lw_q_B <- nb2listw(nb_q, style="B", zero.policy = TRUE)
unlist(spweights.constants(lw_q_B, zero.policy = TRUE))
```

```
##      n     n1     n2     n3     nn     S0     S1     S2 
##    349    348    347    346 121801   1810   3620  41360
```

The Moran's I test results for Massachusetts density between the years 2010-2018 using spatial weights from contiguous neighbors all have a p-value < 2.2e-16 and The Moran Statistic is ~0.62.  


```r
#moran I test for contiguous neighbours
mass_moranI <- list()
for (i in 3:11) {
    mass_moranI[[i-2]] <- capture.output(moran.test(as.numeric(mass_density[[i]]), lw_q_B, zero.policy=TRUE, randomisation = FALSE))
}

for (i in (1:length(mass_moranI))) {print(paste(years[i],": ",mass_moranI[[i]][8], mass_moranI[[i]][9], mass_moranI[[i]][12],mass_moranI[[i]][13])) }
```

```
## [1] "2010-01-01 :  Moran I statistic standard deviate = 19.081, p-value < 2.2e-16 alternative hypothesis: greater       0.625694489      -0.002873563       0.001085180  "
## [1] "2011-01-01 :  Moran I statistic standard deviate = 19.087, p-value < 2.2e-16 alternative hypothesis: greater       0.625888771      -0.002873563       0.001085180  "
## [1] "2012-01-01 :  Moran I statistic standard deviate = 19.116, p-value < 2.2e-16 alternative hypothesis: greater       0.626836683      -0.002873563       0.001085180  "
## [1] "2013-01-01 :  Moran I statistic standard deviate = 19.11, p-value < 2.2e-16 alternative hypothesis: greater       0.626647618      -0.002873563       0.001085180  "
## [1] "2014-01-01 :  Moran I statistic standard deviate = 19.098, p-value < 2.2e-16 alternative hypothesis: greater       0.626270456      -0.002873563       0.001085180  "
## [1] "2015-01-01 :  Moran I statistic standard deviate = 19.067, p-value < 2.2e-16 alternative hypothesis: greater       0.625225455      -0.002873563       0.001085180  "
## [1] "2016-01-01 :  Moran I statistic standard deviate = 19.068, p-value < 2.2e-16 alternative hypothesis: greater       0.625256016      -0.002873563       0.001085180  "
## [1] "2017-01-01 :  Moran I statistic standard deviate = 19.094, p-value < 2.2e-16 alternative hypothesis: greater       0.626127408      -0.002873563       0.001085180  "
## [1] "2018-01-01 :  Moran I statistic standard deviate = 19.133, p-value < 2.2e-16 alternative hypothesis: greater       0.627415241      -0.002873563       0.001085180  "
```

#### Spatial autocorrelation for k-nearest neighbours - Massachusetts Density

The nearest neighbours are calculated from the polygon centroids to ${k}=\sqrt{351}$. Next, the spatial weights are calculated using zero.policy to disclude Massachusetts islands.  


```r
#create k-nearest neighbours
#nearest neighbour using k = sqrt(N)
coords <- st_centroid(st_geometry(mass_density), of_largest_polygon=TRUE)
knn_k <- knearneigh(coords, k=sqrt(351))
nb_k <- knn2nb(knn_k, sym=TRUE)
n_comp <- n.comp.nb(nb_k)
n_comp$nc
```

```
## [1] 1
```

```r
#spatial weights for k-nearest neighbours
lw_q_B_d <- nb2listw(nb_k, style="B", zero.policy = TRUE)
unlist(spweights.constants(lw_q_B_d, zero.policy = TRUE))
```

```
##      n     n1     n2     n3     nn     S0     S1     S2 
##    351    350    349    348 123201   7340  14680 623928
```

The Moran's I test results for Massachusetts density between the years 2010-2018 using spatial weights from k-nearest neighbors all have a p-value < 2.2e-16 and The Moran Statistic is ~0.56.  


```r
#moran I test for k-nearest neighbours
mass_moranI_d <- list()
for (i in 3:11) {
    mass_moranI_d[[i-2]] <- capture.output(moran.test(as.numeric(mass_density[[i]]), lw_q_B_d, zero.policy=TRUE, randomisation = FALSE))
}

for (i in (1:length(mass_moranI_d))) {print(paste(years[i],": ",mass_moranI_d[[i]][7], mass_moranI_d[[i]][8], mass_moranI_d[[i]][11],mass_moranI_d[[i]][12])) }
```

```
## [1] "2010-01-01 :  Moran I statistic standard deviate = 35.613, p-value < 2.2e-16 alternative hypothesis: greater      0.5665928414     -0.0028571429      0.0002556749  "
## [1] "2011-01-01 :  Moran I statistic standard deviate = 35.615, p-value < 2.2e-16 alternative hypothesis: greater      0.5666217699     -0.0028571429      0.0002556749  "
## [1] "2012-01-01 :  Moran I statistic standard deviate = 35.608, p-value < 2.2e-16 alternative hypothesis: greater      0.5665044700     -0.0028571429      0.0002556749  "
## [1] "2013-01-01 :  Moran I statistic standard deviate = 35.539, p-value < 2.2e-16 alternative hypothesis: greater      0.5654020057     -0.0028571429      0.0002556749  "
## [1] "2014-01-01 :  Moran I statistic standard deviate = 35.52, p-value < 2.2e-16 alternative hypothesis: greater      0.5650989112     -0.0028571429      0.0002556749  "
## [1] "2015-01-01 :  Moran I statistic standard deviate = 35.462, p-value < 2.2e-16 alternative hypothesis: greater      0.5641797635     -0.0028571429      0.0002556749  "
## [1] "2016-01-01 :  Moran I statistic standard deviate = 35.412, p-value < 2.2e-16 alternative hypothesis: greater      0.5633747217     -0.0028571429      0.0002556749  "
## [1] "2017-01-01 :  Moran I statistic standard deviate = 35.413, p-value < 2.2e-16 alternative hypothesis: greater      0.5633858582     -0.0028571429      0.0002556749  "
## [1] "2018-01-01 :  Moran I statistic standard deviate = 35.405, p-value < 2.2e-16 alternative hypothesis: greater      0.5632613443     -0.0028571429      0.0002556749  "
```

## Population Density Trends in Most Populated Cities and their Neighboring Towns and Town-Cities

Most Populated Cities in Massachusetts:

1. Boston
2. Worcester
3. Springfield
4. Lowell
5. Cambridge


```r
#filter main cities
m_cities = mass_pop2 %>% filter(Name == "Boston" | Name == "Worcester" | Name == "Springfield"| Name == "Lowell" | Name == "Cambridge")

m_cities_ts = cities_density %>% select("Boston", "Worcester", "Springfield", "Lowell", "Cambridge")
m_cities_ts = ts(m_cities_ts, start = 2010, end = 2018, frequency = 1)
m_cities_ts = cbind(years, m_cities_ts)
m_cities = m_cities%>%select(-geometry,everything())
```

Examinig the population density time series for main series, all cities display a significant increase in population density between 2010-2018. Springfield presents a slight decline between 2014 and 2016 and then the population density continues to increase.


```r
#plot main cities density times series
par(mfrow=c(2,3))
ts.plot(m_cities_ts[,2], main = "Boston Population Density", col = "green", ylab = "Population Density 1/km^2")
ts.plot(m_cities_ts[,3], main = "Worcester Population Density", col = "red", ylab = "Population Density 1/km^2")
ts.plot(m_cities_ts[,4], main = "Springfield Population Density", col = "blue", ylab = "Population Density 1/km^2")
ts.plot(m_cities_ts[,5], main = "Lowell Population Density", col = "orange", ylab = "Population Density 1/km^2")
ts.plot(m_cities_ts[,6], main = "Cambridge Population Density", col = "purple", ylab = "Population Density 1/km^2")
```

![](FinalAssignment_files/figure-html/plot main cities-1.png)<!-- -->

Main cities neighbours plot shows that Boston and Cambridge, two of the main cities in Massachusetts are neighbours of one another.

The density map for main city neighbours presents mainly an increase in density in Boston and Cambridge neighbours.


```r
#find city neighbors - due to an error using r st_intersect the data was collected in QGIS
m_neighbours = read_sf("data/m_neighbours_final.shp", stringsAsFactors = FALSE)
col_names = names(m_cities)
names(m_neighbours) = col_names

par(mfrow=c(1,1))
plot(st_geometry(m_neighbours), col = 'lightblue', main = "Most Populated Cities with their Adjacent Neighbours")
plot(st_geometry(m_cities), col = "red", add = TRUE)
text(st_coordinates(st_centroid(m_cities)), m_cities$Name)
```

![](FinalAssignment_files/figure-html/citie neighbours-1.png)<!-- -->

```r
#plot densities in neighbour cities
plot(m_neighbours[,33:41], key.pos = 1, key.width = lcm(1))
```

![](FinalAssignment_files/figure-html/citie neighbours-2.png)<!-- -->

* Boston Neighbours Population Density Time Series:
The time series present that most of Boston's nieghbours do not have a significant increase in population density over time. Three out of 15 neighbours that do display a slightly more significant increase are Cambridge, Revere and Somerville.

* Worcester Neighbours Population Density Time Series:
The time series present that most of Worcester's nieghbours do not have a significant increase in population density over time. Shrewsbury displays the most significant increase.

* Lowell Neighbours Population Density Time Series:
The time series present that most of Lowell's nieghbours have a slight increase in population density over time. Billerica displays the highest increase.

* Springfield Neighbours Population Density Time Series:
The time series present that most of Springfield's nieghbours do not have an increase in population density over time. Some of the neighbours even show a decline in population density such as, Chicopee.

* Cambridge Neighbours Population Density Time Series:
The time series present that most of Cambridge's nieghbours have an increase in population density over time. The highest increase is presented in Boston, Sommerville and Everett.

![](FinalAssignment_files/figure-html/neighbors time series-1.png)<!-- -->![](FinalAssignment_files/figure-html/neighbors time series-2.png)<!-- -->![](FinalAssignment_files/figure-html/neighbors time series-3.png)<!-- -->![](FinalAssignment_files/figure-html/neighbors time series-4.png)<!-- -->![](FinalAssignment_files/figure-html/neighbors time series-5.png)<!-- -->

#### Spatial autocorrelation for k-nearest neighbours - Main Cities and Their Neighbours Density

The nearest neighbours are calculated from the polygon centroids to ${k}=\sqrt{n}$. The plot of each city and it's neighbours displays the ${k}$-nearest neighbours calculation. Next, the spatial weights are calculated. 



#### Boston and Boston Neighbours spatial autocorrelation

${k}=\sqrt{16}$

Moran's I test results show that in 2010, 2011 and 2012 the Moran I statistic for population density in Boston and Boston neighbours was ~0.3 with p-value = ~0.023. In 2013, 2014, 2015 and 2016 the Moran I Statistic was ~0.29 with p-value = ~0.17. In 2017 and 2018 the Moran Static decreased to ~0.28 with p-value = ~0.017.


```r
#create k-nearest neighbours boston
#k-nearest neighbour k = sqrt(N) = 16
coords_b <- st_centroid(st_geometry(boston), of_largest_polygon=TRUE)
knn_k_b <- knearneigh(coords_b, k=sqrt(16))
nb_k_b <- knn2nb(knn_k_b, sym=TRUE)
n_comp_b <- n.comp.nb(nb_k_b)
n_comp_b$nc
```

```
## [1] 1
```

```r
#plot boston neighbours
plot(st_geometry(boston), border = 'grey', main = "Boston k-nearest Neighbours")
plot(nb_k_b, st_centroid(st_geometry(boston)), pch = 3, cex = .2, add = TRUE, randomisation = FALSE)
```

![](FinalAssignment_files/figure-html/spatial boston, spatial autocorrelation main cities-1.png)<!-- -->

```r
#spatial weights for k-nearest neighbours boston
lw_q_B_d_b <- nb2listw(nb_k_b, style="B", zero.policy = TRUE)
unlist(spweights.constants(lw_q_B_d_b, zero.policy = TRUE))
```

```
##    n   n1   n2   n3   nn   S0   S1   S2 
##   16   15   14   13  256   76  152 1504
```

```r
#moran I test for k-nearest neighbours boston
mass_moranI_d_b <- list()
for (i in 33:41) {
    mass_moranI_d_b[[i-32]] <- capture.output(moran.test(as.numeric(boston[[i]]), lw_q_B_d_b, zero.policy=TRUE, randomisation = FALSE))
}

for (i in (1:length(mass_moranI_d_b))) {print(paste(years[i],": ",mass_moranI_d_b[[i]][7], mass_moranI_d_b[[i]][8], mass_moranI_d_b[[i]][11],mass_moranI_d_b[[i]][12])) }
```

```
## [1] "2010-01-01 :  Moran I statistic standard deviate = 2.829, p-value = 0.002335 alternative hypothesis: greater        0.30651851       -0.06666667        0.01740119  "
## [1] "2011-01-01 :  Moran I statistic standard deviate = 2.8314, p-value = 0.002317 alternative hypothesis: greater        0.30683494       -0.06666667        0.01740119  "
## [1] "2012-01-01 :  Moran I statistic standard deviate = 2.8092, p-value = 0.002483 alternative hypothesis: greater        0.30390078       -0.06666667        0.01740119  "
## [1] "2013-01-01 :  Moran I statistic standard deviate = 2.7624, p-value = 0.002868 alternative hypothesis: greater        0.29773768       -0.06666667        0.01740119  "
## [1] "2014-01-01 :  Moran I statistic standard deviate = 2.7487, p-value = 0.002992 alternative hypothesis: greater        0.29591940       -0.06666667        0.01740119  "
## [1] "2015-01-01 :  Moran I statistic standard deviate = 2.742, p-value = 0.003053 alternative hypothesis: greater        0.29503778       -0.06666667        0.01740119  "
## [1] "2016-01-01 :  Moran I statistic standard deviate = 2.7064, p-value = 0.003401 alternative hypothesis: greater        0.29034191       -0.06666667        0.01740119  "
## [1] "2017-01-01 :  Moran I statistic standard deviate = 2.6846, p-value = 0.003631 alternative hypothesis: greater        0.28746745       -0.06666667        0.01740119  "
## [1] "2018-01-01 :  Moran I statistic standard deviate = 2.6544, p-value = 0.003972 alternative hypothesis: greater        0.28349120       -0.06666667        0.01740119  "
```

#### Worcester and Worcester Neighbours spatial autocorrelation

${k}=\sqrt{9}$

Moran's I test results show that in all the years the Moran I statistic for population density in Worcester and Worcester neighbours was ~-0.39 with p-value = ~0.95.


```r
#create k-nearest neighbours worcester
#k-nearest neighbour k = sqrt(N) = 9
coords_w <- st_centroid(st_geometry(worcester), of_largest_polygon=TRUE)
knn_k_w <- knearneigh(coords_w, k=sqrt(9))
nb_k_w <- knn2nb(knn_k_w, sym=TRUE)
n_comp_w <- n.comp.nb(nb_k_w)
n_comp_w$nc
```

```
## [1] 1
```

```r
#plot worcester neighbours
plot(st_geometry(worcester), border = 'grey', main = "Worcester k-nearest Neighbours")
plot(nb_k_w, st_centroid(st_geometry(worcester)), pch = 3, cex = .2, add = TRUE, randomisation = FALSE)
```

![](FinalAssignment_files/figure-html/spatial worcester, spatial autocorrelation main cities-1.png)<!-- -->

```r
#spatial weights for k-nearest neighbours worcester
lw_q_B_d_w <- nb2listw(nb_k_w, style="B", zero.policy = TRUE)
unlist(spweights.constants(lw_q_B_d_w, zero.policy = TRUE))
```

```
##   n  n1  n2  n3  nn  S0  S1  S2 
##   9   8   7   6  81  32  64 544
```

```r
#moran I test for k-nearest neighbours worcester
mass_moranI_d_w <- list()
for (i in 33:41) {
    mass_moranI_d_w[[i-32]] <- capture.output(moran.test(as.numeric(worcester[[i]]), lw_q_B_d_w, zero.policy=TRUE, randomisation = FALSE))
}

for (i in (1:length(mass_moranI_d_w))) {print(paste(years[i],": ",mass_moranI_d_w[[i]][7], mass_moranI_d_w[[i]][8], mass_moranI_d_w[[i]][11],mass_moranI_d_w[[i]][12])) }
```

```
## [1] "2010-01-01 :  Moran I statistic standard deviate = -1.7182, p-value = 0.9571 alternative hypothesis: greater       -0.39878759       -0.12500000        0.02539062  "
## [1] "2011-01-01 :  Moran I statistic standard deviate = -1.7175, p-value = 0.9571 alternative hypothesis: greater       -0.39867045       -0.12500000        0.02539062  "
## [1] "2012-01-01 :  Moran I statistic standard deviate = -1.714, p-value = 0.9567 alternative hypothesis: greater       -0.39811097       -0.12500000        0.02539062  "
## [1] "2013-01-01 :  Moran I statistic standard deviate = -1.7142, p-value = 0.9567 alternative hypothesis: greater       -0.39814021       -0.12500000        0.02539062  "
## [1] "2014-01-01 :  Moran I statistic standard deviate = -1.7152, p-value = 0.9568 alternative hypothesis: greater       -0.39830964       -0.12500000        0.02539062  "
## [1] "2015-01-01 :  Moran I statistic standard deviate = -1.7135, p-value = 0.9567 alternative hypothesis: greater       -0.39803428       -0.12500000        0.02539062  "
## [1] "2016-01-01 :  Moran I statistic standard deviate = -1.714, p-value = 0.9567 alternative hypothesis: greater       -0.39811668       -0.12500000        0.02539062  "
## [1] "2017-01-01 :  Moran I statistic standard deviate = -1.7059, p-value = 0.956 alternative hypothesis: greater       -0.39682010       -0.12500000        0.02539062  "
## [1] "2018-01-01 :  Moran I statistic standard deviate = -1.6976, p-value = 0.9552 alternative hypothesis: greater       -0.39549689       -0.12500000        0.02539062  "
```

#### Lowell and Lowell Neighbours spatial autocorrelation

${k}=\sqrt{6}$

Moran's I test results show that in 2010, 2011, 2012, 2013, 2014, 2015 and 2016 the Moran I statistic for population density in Lowell and Lowell neighbours was ~-0.44 with p-value = ~0.85. In 2017 and 2018  Moran I statistic was ~-0.44 with p-value = ~0.85.


```r
#create k-nearest neighbours lowell
#k-nearest neighbour k = sqrt(N) = 6
coords_l <- st_centroid(st_geometry(lowell), of_largest_polygon=TRUE)
knn_k_l <- knearneigh(coords_l, k=sqrt(6))
nb_k_l <- knn2nb(knn_k_l, sym=TRUE)
n_comp_l <- n.comp.nb(nb_k_l)
n_comp_l$nc
```

```
## [1] 1
```

```r
#plot lowell neighbours
plot(st_geometry(lowell), border = 'grey', main = "Lowell k-nearest Neighbours")
plot(nb_k_l, st_centroid(st_geometry(lowell)), pch = 3, cex = .2, add = TRUE, randomisation = FALSE)
```

![](FinalAssignment_files/figure-html/spatial lowell, spatial autocorrelation main cities-1.png)<!-- -->

```r
#spatial weights for k-nearest neighbours lowell
lw_q_B_d_l <- nb2listw(nb_k_l, style="B", zero.policy = TRUE)
unlist(spweights.constants(lw_q_B_d_l, zero.policy = TRUE))
```

```
##   n  n1  n2  n3  nn  S0  S1  S2 
##   6   5   4   3  36  16  32 184
```

```r
#moran I test for k-nearest neighbours lowell
mass_moranI_d_l <- list()
for (i in 33:41) {
    mass_moranI_d_l[[i-32]] <- capture.output(moran.test(as.numeric(lowell[[i]]), lw_q_B_d_l, zero.policy=TRUE, randomisation = FALSE))
}

for (i in (1:length(mass_moranI_d_l))) {print(paste(years[i],": ",mass_moranI_d_l[[i]][7], mass_moranI_d_l[[i]][8], mass_moranI_d_l[[i]][11],mass_moranI_d_l[[i]][12])) }
```

```
## [1] "2010-01-01 :  Moran I statistic standard deviate = -1.0546, p-value = 0.8542 alternative hypothesis: greater       -0.43833430       -0.20000000        0.05107143  "
## [1] "2011-01-01 :  Moran I statistic standard deviate = -1.0569, p-value = 0.8547 alternative hypothesis: greater       -0.43883854       -0.20000000        0.05107143  "
## [1] "2012-01-01 :  Moran I statistic standard deviate = -1.0554, p-value = 0.8544 alternative hypothesis: greater       -0.43851108       -0.20000000        0.05107143  "
## [1] "2013-01-01 :  Moran I statistic standard deviate = -1.0581, p-value = 0.855 alternative hypothesis: greater       -0.43912302       -0.20000000        0.05107143  "
## [1] "2014-01-01 :  Moran I statistic standard deviate = -1.0583, p-value = 0.855 alternative hypothesis: greater       -0.43915791       -0.20000000        0.05107143  "
## [1] "2015-01-01 :  Moran I statistic standard deviate = -1.0576, p-value = 0.8549 alternative hypothesis: greater       -0.43901681       -0.20000000        0.05107143  "
## [1] "2016-01-01 :  Moran I statistic standard deviate = -1.0574, p-value = 0.8548 alternative hypothesis: greater       -0.43895943       -0.20000000        0.05107143  "
## [1] "2017-01-01 :  Moran I statistic standard deviate = -1.0671, p-value = 0.857 alternative hypothesis: greater       -0.44115926       -0.20000000        0.05107143  "
## [1] "2018-01-01 :  Moran I statistic standard deviate = -1.0652, p-value = 0.8566 alternative hypothesis: greater       -0.44073194       -0.20000000        0.05107143  "
```

#### Springfield and Springfield Neighbours spatial autocorrelation

${k}=\sqrt{8}$

Moran's I test results show that in all the years the Moran I statistic for population density in Springfield and Springfield neighbours was ~-0.18 with p-value = ~0.56. 


```r
#create k-nearest neighbours springfield
#k-nearest neighbour k = sqrt(N) = 8
coords_s <- st_centroid(st_geometry(springfield), of_largest_polygon=TRUE)
knn_k_s <- knearneigh(coords_s, k=sqrt(8))
nb_k_s <- knn2nb(knn_k_s, sym=TRUE)
n_comp_s <- n.comp.nb(nb_k_s)
n_comp_s$nc
```

```
## [1] 1
```

```r
#plot springfield neighbours
plot(st_geometry(springfield), border = 'grey', main = "Springfield k-nearest Neighbours")
plot(nb_k_s, st_centroid(st_geometry(springfield)), pch = 3, cex = .2, add = TRUE, randomisation = FALSE)
```

![](FinalAssignment_files/figure-html/spatial springfield, spatial autocorrelation main cities-1.png)<!-- -->

```r
#spatial weights for k-nearest neighbours springfield
lw_q_B_d_s <- nb2listw(nb_k_s, style="B", zero.policy = TRUE)
unlist(spweights.constants(lw_q_B_d_s, zero.policy = TRUE))
```

```
##   n  n1  n2  n3  nn  S0  S1  S2 
##   8   7   6   5  64  18  36 168
```

```r
#moran I test for k-nearest neighbours springfield
mass_moranI_d_s <- list()
for (i in 33:41) {
    mass_moranI_d_s[[i-32]] <- capture.output(moran.test(as.numeric(springfield[[i]]), lw_q_B_d_s, zero.policy=TRUE, randomisation = FALSE))
}

for (i in (1:length(mass_moranI_d_s))) {print(paste(years[i],": ",mass_moranI_d_s[[i]][7], mass_moranI_d_s[[i]][8], mass_moranI_d_s[[i]][11],mass_moranI_d_s[[i]][12])) }
```

```
## [1] "2010-01-01 :  Moran I statistic standard deviate = -0.17193, p-value = 0.5683 alternative hypothesis: greater       -0.18970270       -0.14285714        0.07424204  "
## [1] "2011-01-01 :  Moran I statistic standard deviate = -0.16823, p-value = 0.5668 alternative hypothesis: greater       -0.18869574       -0.14285714        0.07424204  "
## [1] "2012-01-01 :  Moran I statistic standard deviate = -0.16204, p-value = 0.5644 alternative hypothesis: greater       -0.18700987       -0.14285714        0.07424204  "
## [1] "2013-01-01 :  Moran I statistic standard deviate = -0.15713, p-value = 0.5624 alternative hypothesis: greater       -0.18566989       -0.14285714        0.07424204  "
## [1] "2014-01-01 :  Moran I statistic standard deviate = -0.14723, p-value = 0.5585 alternative hypothesis: greater       -0.18297452       -0.14285714        0.07424204  "
## [1] "2015-01-01 :  Moran I statistic standard deviate = -0.15321, p-value = 0.5609 alternative hypothesis: greater       -0.18460178       -0.14285714        0.07424204  "
## [1] "2016-01-01 :  Moran I statistic standard deviate = -0.157, p-value = 0.5624 alternative hypothesis: greater       -0.18563538       -0.14285714        0.07424204  "
## [1] "2017-01-01 :  Moran I statistic standard deviate = -0.15144, p-value = 0.5602 alternative hypothesis: greater       -0.18412022       -0.14285714        0.07424204  "
## [1] "2018-01-01 :  Moran I statistic standard deviate = -0.1523, p-value = 0.5605 alternative hypothesis: greater       -0.18435357       -0.14285714        0.07424204  "
```

#### Cambridge and Cambridge Neighbours spatial autocorrelation

${k}=\sqrt{7}$

Moran's I test results show that in 2010 Moran I statistic for population density in Cambridge and Cambridge neighbours was ~0.089 with p-value = ~0.15. In 2011, the Moran I Statistic was ~0.087 with p-value = ~0.15. In 2012, the Moran I Statistic was ~0.09 with p-value = ~0.15. In 2013, the Moran I Statistic was ~0.096 with p-value = ~0.14. In 2014, the Moran I Statistic was ~0.1 with p-value = ~0.14. In 2015, the Moran I Statistic was ~0.099 with p-value = ~0.13. In 2016, the Moran I Statistic was ~0.1 with p-value = ~0.13. In 2017, the Moran I Statistic was ~0.105 with p-value = ~0.13, and in 2018 the Moran Static was ~0.107 with p-value = ~0.134.


```r
#create k-nearest neighbours cambridge
#k-nearest neighbour k = sqrt(N) = 7
coords_c <- st_centroid(st_geometry(cambridge), of_largest_polygon=TRUE)
knn_k_c <- knearneigh(coords_c, k=sqrt(7))
nb_k_c <- knn2nb(knn_k_c, sym=TRUE)
n_comp_c <- n.comp.nb(nb_k_c)
n_comp_c$nc
```

```
## [1] 1
```

```r
#plot cambridge neighbours
plot(st_geometry(cambridge), border = 'grey', main = "Cambridge k-nearest Neighbours")
plot(nb_k_c, st_centroid(st_geometry(cambridge)), pch = 3, cex = .2, add = TRUE, randomisation = FALSE)
```

![](FinalAssignment_files/figure-html/spatial cambridge, spatial autocorrelation main cities-1.png)<!-- -->

```r
#spatial weights for k-nearest neighbours cambridge
lw_q_B_d_c <- nb2listw(nb_k_c, style="B", zero.policy = TRUE)
unlist(spweights.constants(lw_q_B_d_c, zero.policy = TRUE))
```

```
##   n  n1  n2  n3  nn  S0  S1  S2 
##   7   6   5   4  49  18  36 192
```

```r
#moran I test for k-nearest neighbours cambridge
mass_moranI_d_c <- list()
for (i in 33:41) {
    mass_moranI_d_c[[i-32]] <- capture.output(moran.test(as.numeric(cambridge[[i]]), lw_q_B_d_c, zero.policy=TRUE, randomisation = FALSE))
}

for (i in (1:length(mass_moranI_d_c))) {print(paste(years[i],": ",mass_moranI_d_c[[i]][7], mass_moranI_d_c[[i]][8], mass_moranI_d_c[[i]][11],mass_moranI_d_c[[i]][12])) }
```

```
## [1] "2010-01-01 :  Moran I statistic standard deviate = 1.029, p-value = 0.1517 alternative hypothesis: greater        0.08900146       -0.16666667        0.06172840  "
## [1] "2011-01-01 :  Moran I statistic standard deviate = 1.023, p-value = 0.1532 alternative hypothesis: greater        0.08749731       -0.16666667        0.06172840  "
## [1] "2012-01-01 :  Moran I statistic standard deviate = 1.0341, p-value = 0.1505 alternative hypothesis: greater        0.09025423       -0.16666667        0.06172840  "
## [1] "2013-01-01 :  Moran I statistic standard deviate = 1.0583, p-value = 0.145 alternative hypothesis: greater        0.09627147       -0.16666667        0.06172840  "
## [1] "2014-01-01 :  Moran I statistic standard deviate = 1.0755, p-value = 0.1411 alternative hypothesis: greater         0.1005463        -0.1666667         0.0617284  "
## [1] "2015-01-01 :  Moran I statistic standard deviate = 1.0716, p-value = 0.1419 alternative hypothesis: greater        0.09958097       -0.16666667        0.06172840  "
## [1] "2016-01-01 :  Moran I statistic standard deviate = 1.0802, p-value = 0.14 alternative hypothesis: greater         0.1017117        -0.1666667         0.0617284  "
## [1] "2017-01-01 :  Moran I statistic standard deviate = 1.094, p-value = 0.137 alternative hypothesis: greater         0.1051380        -0.1666667         0.0617284  "
## [1] "2018-01-01 :  Moran I statistic standard deviate = 1.1039, p-value = 0.1348 alternative hypothesis: greater         0.1076056        -0.1666667         0.0617284  "
```

***

## Discussion and Conclusions

The population growth analysis for Massachusetts settlement types presented an unexpected decline in growth rate when analyzing cities. Although the growth rate between 2017 and 2018 is begining to increase again. Population growth density analysis for Massachusetts settlement types displayed a similar result to the population growth analysis. In addition the density population growth rate in towns showed an increase between 2017 and 2018. This may be a result of urban sprawl and increasing population in total. 

The spatial autocorrelation for Massachusetts population density results display a higher significant spatial autocorrelation between contiguous neighbours based spatial weight matrix compared to the ${k}$-nearest based spatial weight matrix. In contrast to the paper about Spatial Regression Models for Demographic Analysis by Chi & Zhu (2007). 
The results of Moran's I test confirm that population density at one place has a significant effect on population density at another place and the concentration of population is not randomly dispersed. Which correspond to the hypothesis that population in the years 2010 to 2018 in Massachusetts is becoming more dense around a city and not just at a random location.

Population density in main cities in Massachusetts is increasing significantly and their neighbours population density is varied. However, most of the neighbours are not presenting a substantial increase in population density. The highest increase is occuring in Boston and Cambridge and their neighbours. 

The spatial autocorrelation in main cities presented only a significant spatial autocorraltion bewteen Boston and it's neighbours. The Moran's I test result was not as high when comparing the results to population density in all of Massachusetts. Nevertheless, there is an effect on population density between Boston and it's neighbours and the population is not randomly dispersed. This result coincides with the hypothesis that population is becoming more dense in Massachusetts cities the years 2010 to 2018.

This study is important for better understanding the changes and trends of population density in urban and suburban areas. 
This study finds that there is a significant spatial autocorrelation in Massachusetts population density and an increase in population density between 2010-2018. In this study the population density change through a small range of years was used, this allowed only to roughly examine the spatio-temporal population density changes. For further analysis a wider temporal range should be evaluated. 


## References
* Bivand, R.(2019).Creating Neighbours.Retrieved from https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
* Chi, G. & Zhu, J.(2007).Spatial Regression Models for Demographic Analysis.Popul Res Policy Rev (2008) 27:17-42.doi:10.1007/s11113-007-9051-8
* Cox, W.(July 16, 2015).The Evolving Urban Form: Sprawling Boston. Retrieved from https://www.newgeography.com/content/004987-the-evolving-urban-form-sprawling-boston
* Federal Research Division of the Library of Congress(n.d).Country Studies.Retrieved from http://countrystudies.us/
* Moran's I: Definition, Examples. Retrieved from https://www.statisticshowto.datasciencecentral.com/morans-i/
* Paige, W.S., Ryan, R.L., Lerman, S.B., & Tooke, K.A.,(2011).Social and institutional factors associated with land use and forest conservation along two urban gradients in Massachusetts.Landscape and Urban Planning (102) 82- 92.doi:dx.doi.org/10.1016/j.landurbplan.2011.03.012
* Pebesma, E. & Bivand, R.(2020).Spatial Data Science.Retrieved from https://keen-swartz-3146c4.netlify.com/index.html
* Subramanian, D.(June 8, 2019).A Simple Introduction to K-Nearest Neighbors Algorithm.Retrieved from https://towardsdatascience.com/a-simple-introduction-to-k-nearest-neighbors-algorithm-b3519ed98e
* United Nations - Department of Economic and Social Affairs Population Dynamics(2019).2018 Revision of World Urbanization Prospects.Retrieved from https://population.un.org/wup/


