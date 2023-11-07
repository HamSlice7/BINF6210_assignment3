#1.Load the necessary libraries----
library("tidyverse")
library("vegan")
library("fossil")
library("ggplot2")
library("ggmap")

#2.Load in the Nymphalidae tsv----
#read in the Nymphalidae data frame from BOLD
df_n <- read_tsv('http://www.boldsystems.org/index.php/API_Public/combined?taxon=Nymphalidae&format=tsv')
#3.Looking at the number of samples at each latitudinal range----

#The function serves as a tool for exploratory analysis of my data. It selects sampleid and lat columns from df_n and filters the desired latitudinal range. It then groups the data by latitude and counts the number of samples at each latitude group. It then pulls the count column and then calculates the sum of this column to get the number of samples made of the taxa at the inputed latitudinal band
num_of_samples_nymph <- function(minimum_latitude, max_latitude) {
  sample_count <- df_n %>%
    select(sampleid, lat) %>%
    filter(!is.na(lat) & lat >= minimum_latitude & lat < max_latitude) %>%
    group_by(lat) %>%
    summarise(count = n()) %>%
    pull(count) %>%
    sum() 
  return(sample_count) }

#creating a vector of function calls to be used in the nymph_sampling data frame.
nymph_sampling_effort <- c(num_of_samples_nymph(0, 15), num_of_samples_nymph(16, 30), num_of_samples_nymph(31, 45), num_of_samples_nymph(46, 60), num_of_samples_nymph(61, 75), num_of_samples_nymph(76, 90))

#creating a vector containing my latitudinal bands to be used in the nymph_sampling data frame.
latitude <- c('0-15', '16-30', '31-45', '46-60', '61-75', '76-90')

#creating nymph_sampling data frame from latitude and nymph_sampling_effort
nymph_sampling <- data.frame(latitude, nymph_sampling_effort)
#4.Creating a map of sampling efforts----


#Sample data frame with latitude and longitude 
df_n_map <- df_n %>%
  select(sampleid, lat, lon) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  filter(lat >= 0)

register_google(key = "AIzaSyBxWwRVQigeDNM_MhkZlR2fhXHN8e8Vqbk")

#create map plot
map_nymph <- ggmap(get_googlemap(center = c(lon = 0, lat = 40), zoom = 1, maptype = 'terrain', color = 'color')) +
  geom_point(data = df_n_map, aes(x=lon, y=lat), color = 'red', size = 0.1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
  ) 

print(map_nymph)
#5.Creating a function to calculate species richness using taxonomic id and vegan----

#Here, I create a function that will select species_name and lat and filter out any na's. It will then filter the desired latitudinal range of the user. Next, it will then group by species_name and count the number of species_name groups. It will then reformat the data so that it can be used for the specnumber() function in the vegan package. The function will then return the estimated species richness.
speciesrichness_taxid_vegan_nymph <- function(min_lat, max_lat) {
  c_n <- df_n %>%
    select(species_name, lat) %>%
    filter(!is.na(species_name), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(species_name) %>%
    count(species_name) %>%
    pivot_wider(names_from = species_name, values_from = n) %>%
    vegan::specnumber() %>%
    return(c_n) 
}

#calling the function based on my desired latitudinal ranges
n_taxid_vegan_015 <- speciesrichness_taxid_vegan_nymph(0, 15)
n_taxid_vegan_1630 <- speciesrichness_taxid_vegan_nymph(16, 30)
n_taxid_vegan_3145 <- speciesrichness_taxid_vegan_nymph(31, 45)
n_taxid_vegan_4660 <- speciesrichness_taxid_vegan_nymph(46, 60)
n_taxid_vegan_6175 <- speciesrichness_taxid_vegan_nymph(61, 75)
n_taxid_vegan_7690 <- speciesrichness_taxid_vegan_nymph(76, 90)




#6.Creating a function to calculate species richness using bin id and vegan----

#Here, I create a function that will select bin_uri and lat and filter out any NA's. It will then filter the desired latitudinal range of the user. Next, it will then group by bin_uri and count the number of bin_uri groups. It will then reformat the data so that it can be used for the specnumber() function in the vegan package. The function will then return the estimated species richness.
speciesrichness_binid_vegan_nymph <- function(min_lat, max_lat) {
  c_n <- df_n %>%
    select(bin_uri, lat) %>%
    filter(!is.na(bin_uri), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(bin_uri) %>%
    count(bin_uri) %>%
    pivot_wider(names_from = bin_uri, values_from = n) %>%
    vegan::specnumber() %>%
    return(c_n) 
}

#calling the function based on my desired latitudinal ranges
n_binid_vegan_015 <- speciesrichness_binid_vegan_nymph(0, 15)
n_binid_vegan_1630 <- speciesrichness_binid_vegan_nymph(16, 30)
n_binid_vegan_3145 <- speciesrichness_binid_vegan_nymph(31, 45)
n_binid_vegan_4660 <- speciesrichness_binid_vegan_nymph(46, 60)
n_binid_vegan_6175 <- speciesrichness_binid_vegan_nymph(61, 75)
n_binid_vegan_7690 <- speciesrichness_binid_vegan_nymph(76, 90)



#7.Creating a function to calculate species richness using taxonomic id and chao1----

#Here, I create a function that will select species_name and lat and filter out any na's. It will then filter the desired latitudinal range of the user. Next, it will then group by species_name and count the number of species_name groups. It will then reformat the data so that it can be used for the chao1 function in the fossil package. The function will then return the estimated species richness.
speciesrichness_taxid_chao1_nymph <- function(min_lat, max_lat) {
  c_n <- df_n %>%
    select(species_name, lat) %>%
    filter(!is.na(species_name), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(species_name) %>%
    count(species_name) %>%
    pivot_wider(names_from = species_name, values_from = n) %>%
    fossil::chao1(taxa.row = TRUE) %>%
    return(c_n) 
}


#calling function with latitdudinal bands as inputs
n_taxid_chao1_015 <- speciesrichness_taxid_chao1_nymph(0, 15)
n_taxid_chao1_1630 <- speciesrichness_taxid_chao1_nymph(16, 30)
n_taxid_chao1_3145 <- speciesrichness_taxid_chao1_nymph(31, 45)
n_taxid_chao1_4660 <- speciesrichness_taxid_chao1_nymph(46, 60)
n_taxid_chao1_6175 <- speciesrichness_taxid_chao1_nymph(61, 75)
n_taxid_chao1_7690 <- speciesrichness_taxid_chao1_nymph(76, 90)

#8.Creating a function to calculate species richness using bin id and chao1----

#Here, I create a function that will select bin_uri and lat and filter out any NA's. It will then filter the desired latitudinal range of the user. Next, it will then group by bin_uri and count the number of bin_uri groups. It will then reformat the data so that it can be used for the chao1() function in the fossil package. The function will then return the estimated species richness.
speciesrichness_binid_chao1_nymph <- function(min_lat, max_lat) {
  c_n <- df_n %>%
    select(bin_uri, lat) %>%
    filter(!is.na(bin_uri), !is.na(lat)) %>%
    filter(lat >= min_lat & lat < max_lat) %>%
    group_by(bin_uri) %>%
    count(bin_uri) %>%
    pivot_wider(names_from = bin_uri, values_from = n) %>%
    fossil::chao1(taxa.row = TRUE) %>%
    return(c_n) 
}

#calling function with latitudinal bands as inputs
n_binid_chao1_015 <- speciesrichness_binid_chao1_nymph(0, 15)
n_binid_chao1_1630 <- speciesrichness_binid_chao1_nymph(16, 30)
n_binid_chao1_3145 <- speciesrichness_binid_chao1_nymph(31, 45)
n_binid_chao1_4660 <- speciesrichness_binid_chao1_nymph(46, 60)
n_binid_chao1_6175 <- speciesrichness_binid_chao1_nymph(61, 75)
n_binid_chao1_7690 <- speciesrichness_binid_chao1_nymph(76, 90)




#9.Creating a group bar plot to visualize species richness and calculation groups across latitude----

#organizing my data into the proper data frame format to make it compatible with a group bar chart
latitudes <- c('0-15','0-15', '0-15', '0-15', '16-30','16-30', '16-30', '16-30', '31-45','31-45', '31-45', '31-45','46-60','46-60', '46-60', '46-60', '61-75','61-75', '61-75', '61-75', '76-90', '76-90', '76-90', '76-90')

calculation_groups <- c("Chao1  TaxaID", "Vegan TaxaID", "Chao1 BinID", "Vegan BinID")

species_richness_calcgroups_nymph <- c(n_taxid_chao1_015, n_taxid_vegan_015,n_binid_chao1_015, n_binid_vegan_015, n_taxid_chao1_1630, n_taxid_vegan_1630,n_binid_chao1_1630, n_binid_vegan_1630, n_taxid_chao1_3145, n_taxid_vegan_3145,n_binid_chao1_3145, n_binid_vegan_3145, n_taxid_chao1_4660, n_taxid_vegan_4660, n_binid_chao1_4660, n_binid_vegan_4660, n_taxid_chao1_6175,n_taxid_vegan_6175,n_binid_chao1_6175, n_binid_vegan_6175, n_taxid_chao1_7690,n_taxid_vegan_7690, n_binid_chao1_7690, n_binid_vegan_7690)


#creating the data frame
df_nymph_species_richness_calculations <- data.frame(latitudes, calculation_groups, species_richness_calcgroups_nymph)


#creating a bar plot to visualize the results across the calculations. I am using ggplot to create a grouped bar char of my data so that I can visualize each species richness estimate at each of the latitudinal bands. I am ensuring the bar plot colours are colour blind friendly by using 'Dark2' palette.
Nymph_plot <- df_nymph_species_richness_calculations %>%
  ggplot(aes(x = latitudes, y = species_richness_calcgroups_nymph, fill = calculation_groups ))+
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.7)+
  labs(
    title = 'Species Richness Estimates of Nymphalidae and Latitude',
    x = 'Latitudes',
    y = 'Estimated Species Richness',
    fill = 'Calculations') +
  scale_fill_brewer(palette = 'Dark2') +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(size = 9)) +
  coord_flip()

Nymph_plot


















