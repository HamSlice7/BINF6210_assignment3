#1.load the necessary libraries----
library("tidyverse")
library("vegan")
library("fossil")
library("ggmap")
library("ggplot2")


#2. Load in Aves tsv----

#read in the aves data frame from BOLD
df_a <- read_tsv('"http://www.boldsystems.org/index.php/API_Public/combined?taxon=Aves&format=tsv')


#3.Number of samples based on latitude----

#The function serves as a tool for exploratory analysis of my data. It selects sampleid and lat columns from df_a and filters the desired latitudinal range. It then groups the data by latitude and counts the number of samples at each latitude group. It then pulls the count column and then calculates the sum of this column to get the number of samples made of the taxa at the inputed latitudinal band
num_of_samples_aves <- function(minimum_latitude, max_latitude) {
  sample_count <- df_a %>%
    select(sampleid, lat) %>%
    filter(!is.na(lat) & lat >= minimum_latitude & lat < max_latitude) %>%
    group_by(lat) %>%
    summarise(count = n()) %>%
    pull(count) %>%
    sum() %>%
    return(sample_count) }

aves_sampling_effort <- c(num_of_samples_aves(0, 15), num_of_samples_aves(16, 30), num_of_samples_aves(31, 45), num_of_samples_aves(46, 60), num_of_samples_aves(61, 75), num_of_samples_aves(76, 90))

latitude <- c('0-15', '16-30', '31-45', '46-60', '61-75', '76-90')

aves_sampling <- data.frame(latitude, aves_sampling_effort)
#4. Creating a map of sampling efforts for Aves----


#Sample data frame with latitude and longitude 
df_a_map <- df_a %>%
  select(sampleid, lat, lon) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  filter(lat >= 0)

register_google(key = "AIzaSyBxWwRVQigeDNM_MhkZlR2fhXHN8e8Vqbk")

#create map plot of sampling efforts for Aves
map_aves <- ggmap(get_googlemap(center = c(lon = 0, lat = 40), zoom = 1, maptype = 'terrain', color = 'color')) +
  geom_point(data = df_a_map, aes(x=lon, y=lat), color = 'blue', size = 0.1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  ) 

print(map_aves)

#5.Function that estimates species richness using vegan and taxonomic identification   ----

#Here, I create a function that will select species_name and lat and filter out any na's. It will then filter the desired latitudinal range of the user. Next, it will then group by species_name and count the number of species_name groups. It will then reformat the data so that it can be used for the specnumber() function in the vegan package. The function will then return the estimated species richness.
speciesrichness_name_f <- function(min_lat, max_lat) {
  c_a <- df_a %>%
    select(species_name, lat) %>%
    filter(!is.na(species_name), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(species_name) %>%
    count(species_name) %>%
    pivot_wider(names_from = species_name, values_from = n) %>%
    vegan::specnumber() %>%
    return(c_a) 
}
#calling the function with my latitudinal bands as input
name_vegan_015 <- speciesrichness_name_f(0, 15)
name_vegan_1630 <- speciesrichness_name_f(16, 30)
name_vegan_3145 <- speciesrichness_name_f(31, 45)
name_vegan_4660 <- speciesrichness_name_f(46, 60)
name_vegan_6175 <- speciesrichness_name_f(61, 75)
name_vegan_7690 <- speciesrichness_name_f(76, 90)



#6.Function that estimates species richness using chao1 and bin identification----

#Here, I create a function that will select bin_uri and lat and filter out any NA's. It will then filter the desired latitudinal range of the user. Next, it will then group by bin_uri and count the number of bin_uri groups. It will then reformat the data so that it can be used for the chao1 function in the fossil package. The function will then return the estimated species richness.
chao1f <- function(min_lat, max_lat) {
  c_a <- df_a %>%
    select(bin_uri, lat) %>%
    filter(!is.na(bin_uri), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(bin_uri) %>%
    count(bin_uri) %>%
    pivot_wider(names_from = bin_uri, values_from = n) %>%
    fossil::chao1(taxa.row = TRUE) %>%
    return(c_a) 
  }

#calling the function with my latitudinal bands as input
chao1f_015 <- chao1f(0, 15)
chao1f_1630 <- chao1f(16, 30)
chao1f_3145 <- chao1f(31, 45)
chao1f_4660 <- chao1f(46, 60)
chao1f_6175 <- chao1f(61, 75)
chao1f_7690 <- chao1f(76, 90)




    


#7. Function that estimates species richness using vegan and bin identification----

#Here, I create a function that will select bin_uri and lat and filter out any NA's. It will then filter the desired latitudinal range of the user. Next, it will then group by bin_uri and count the number of bin_uri groups. It will then reformat the data so that it can be used for the specnumber() function in the vegan package. The function will then return the estimated species richness.
bin_speciesrichness_f <- function(min_lat, max_lat) {
  c_a <- df_a %>%
    select(bin_uri, lat) %>%
    filter(!is.na(bin_uri), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(bin_uri) %>%
    count(bin_uri) %>%
    pivot_wider(names_from = bin_uri, values_from = n) %>%
    vegan::specnumber() %>%
    return(c_a) 
}
#calling the function with latitudinal bands as inputs
bin_vegan_015 <- bin_speciesrichness_f(0, 15)
bin_vegan_1630 <- bin_speciesrichness_f(16, 30)
bin_vegan_3145 <- bin_speciesrichness_f(31, 45)
bin_vegan_4660 <- bin_speciesrichness_f(46, 60)
bin_vegan_6175 <- bin_speciesrichness_f(61, 75)
bin_vegan_7690 <- bin_speciesrichness_f(76, 90)



#8. Function that estimates species richness using Choa1 and taxanomic identification----

#Here, I create a function that will select species_name and lat and filter out any na's. It will then filter the desired latitudinal range of the user. Next, it will then group by species_name and count the number of species_name groups. It will then reformat the data so that it can be used for the chao1 function in the fossil package. The function will then return the estimated species richness.
chao1_name_f <- function(min_lat, max_lat) {
  c_a <- df_a %>%
    select(species_name, lat) %>%
    filter(!is.na(species_name), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(species_name) %>%
    count(species_name) %>%
    pivot_wider(names_from = species_name, values_from = n) %>%
    fossil::chao1(taxa.row = TRUE) %>%
    return(c_a) 
}

#calling function with latitudinal bands as inputs
chao1_namef_015 <- chao1_name_f(0, 15)
chao1_namef_1630 <- chao1_name_f(16, 30)
chao1_namef_3145 <- chao1_name_f(31, 45)
chao1_namef_4660 <- chao1_name_f(46, 60)
chao1_namef_6175 <- chao1_name_f(61, 75)
chao1_namef_7690 <- chao1_name_f(76, 90)







#9.creating a bar plot to visualize the results across the calculations.----
#creating a data frame to organize my results across the calculations
calculation_groups <- c("Chao1  TaxaID", "Vegan TaxaID", "Chao1 BinID", "Vegan BinID")

species_richness_aves <- c(chao1_namef_015, name_vegan_015,chao1f_015,bin_vegan_015, chao1_namef_1630,name_vegan_1630,chao1f_1630, bin_vegan_1630, chao1_namef_3145,name_vegan_3145,chao1f_3145, bin_vegan_3145, chao1_namef_4660, name_vegan_4660, chao1f_4660,bin_vegan_4660, chao1_namef_6175,name_vegan_6175,bin_vegan_6175, chao1f_6175, chao1_namef_7690,name_vegan_7690, chao1f_7690, bin_vegan_7690)
#creating the data frame
df_aves_species_richness_calculations <- data.frame(latitude, calculation_groups, species_richness_aves)
#creating a bar plot to visualize the results across the calculations. I am using ggplot to create a grouped bar char of my data so that I can visualize each species richness estimate at each of the latitudinal bands. I am ensuring the bar plot colours are colour blind friendly by using 'Dark2' palette.
aves_plot <- df_aves_species_richness_calculations %>%
  ggplot(aes(x = latitudes, y = species_richness_aves, fill = calculation_groups))+
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.7)+
  labs(
    title = 'Species Richness Estimates of Aves and Latitude',
    x = 'Latitudes',
    y = 'Species Richness Estimates',
    fill = 'Calculations') +
  scale_fill_brewer(palette = 'Dark2') +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 9)) +
  coord_flip()

aves_plot
  
