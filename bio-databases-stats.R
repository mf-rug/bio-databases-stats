library(rvest)
library(tidyverse)
library(readxl)
library(RCurl)
library(XML)
library(scales)
library(patchwork)

# Genbank
url <- "https://www.ncbi.nlm.nih.gov/genbank/statistics/"  # Replace with the actual URL
webpage <- read_html(url)

tables <- html_table(webpage, fill = TRUE, header = FALSE)
# Assuming the table you want is the first one
genbank.df <- tables[[1]]

colnames(genbank.df) <- apply(genbank.df, 2, function(x) paste0(x[1], x[2]))
genbank.df <- genbank.df[3:nrow(genbank.df),]

genbank.df.l <- pivot_longer(genbank.df[, c('Date', colnames(genbank.df)[str_detect(colnames(genbank.df), 'Sequences$')])], cols = -1, names_to = 'database', values_to = 'entries')

genbank.df.l[genbank.df.l$entries == '', 'entries'] <- '0'
genbank.df.l$entries <- as.numeric(genbank.df.l$entries)
genbank.df.l$Date <- as.Date(paste("01", genbank.df.l$Date), format="%d %b %Y")

###



#uniprot historic
zip_url <- "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/43/D1/10.1093_nar_gku989/2/gku989_Supplementary_Data.zip?Expires=1711980067&Signature=ojSy741p5oxZYBAEdIgwHPjUrR5b5yi7GRdSVuVDWBvKpfh2Vk9XEGTXsClxh73ClI5kQK8g0DAHSs5in52KAhIsG3lZAKAJFQYO1TAqxXGj~HsvTTFniKnA8oJKAxdGl~aG2TJyXZ8R~e9fGpn8UUgQniQxKICqgqhkg1LOn4HIC~aoxB3~Bid8nk-t1iqBO~yQDSXCrafV~2nh80QS9SBwmh4I3k6-9CbIt7KzeJkfxeZhFUbYhhBH6mPi89BYXQ3Oj46q9g3cDVtSo7TytT0H5G5H3GaxPLydjDctmOtP67seYvPojRja4ahneTIJVYfcs1K4b-HhwykCU0NshA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA" 
zip_file <- tempfile()  # Creates a temporary file
download.file(zip_url, zip_file)

temp_dir <- tempdir()  # Get or create a temporary directory
unzip(zip_file, exdir=temp_dir)

files <- list.files(temp_dir, pattern = "\\.xlsx$", full.names = TRUE)
excel_file <- files[1]  # Assuming it's the first Excel file in the list

uniprot.df <- read_excel(excel_file)
colnames(uniprot.df) <- uniprot.df[26,]
uniprot.df <- uniprot.df[27:205,]
uniprot.df <- uniprot.df %>%
  mutate(Date = as.Date(as.numeric(Date), origin = "1899-12-30"))
# uniprot.df$Date <- as.Date(paste0(uniprot.df$Release, "_01"), format="%Y_%m_%d")
uniprot.df$UniProtKB <- as.numeric(uniprot.df$UniProtKB)

# uniprot.df.l <- pivot_longer(uniprot.df[,-1], col = seq(2, ncol(uniprot.df) -1),names_to = 'database', values_to = 'entries')
# uniprot.df.l$entries <- as.numeric(uniprot.df.l$entries)


# uniprot new
uniprot_new_folders <- getURL('https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/' , dirlistonly = TRUE) %>% getHTMLLinks()
uniprot_new_folders <- uniprot_new_folders[str_detect(uniprot_new_folders, '^release-[0-9]{4}_[0-9]{2}/$')]

uniprot_new.df <- data.frame()
for (i in seq_along(uniprot_new_folders)) {
  upurl1 <- paste0('https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/', uniprot_new_folders[i], '/knowledgebase/UniProtKB_TrEMBL-relstat.html')
  print(upurl1)
  uniprot_recent <- try(read_html(upurl1), silent = TRUE)
  
  # Check if the first attempt was successful
  if (inherits(uniprot_recent, "try-error")) {
    # Second URL to check if the first one fails
    upurl2 <- paste0('https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/', uniprot_new_folders[i], '/knowledgebase/TrEMBL_statistics.html')
    uniprot_recent <- try(read_html(upurl2), silent = TRUE)
    
    # Check if the second attempt was also unsuccessful
    if (inherits(uniprot_recent, "try-error")) {
      next  # Skip to the next iteration if both attempts fail
    }
  }
  uniprot_recent <- uniprot_recent %>% html_text() %>% str_extract('(?<=contains )[0-9]+(?= sequence)') %>% as.numeric()
  if (!is.na(as.numeric(uniprot_recent))) {
    uniprot_new.df[i, 'Date'] <- as.Date(paste0(str_extract(uniprot_new_folders[i], '(?<=release-)[0-9]{4}_[0-9]{2}'), '_01'), format="%Y_%m_%d")
    uniprot_new.df[i, 'UniProtKB'] <- as.numeric(uniprot_recent)
  }
}

uniprot_new.df <- uniprot_new.df[!is.na(uniprot_new.df$Date),]

# merge new and old
uniprot.df <- uniprot.df[!uniprot.df$Date %in% uniprot_new.df$Date,]
uniprot_full.df <- rbind(uniprot.df[, c('Date', 'UniProtKB')], uniprot_new.df)

yearsU <- range(format(uniprot_full.df$Date, "%Y"))

years <- c('1990', max(format(genbank.df.l$Date, "%Y")))
jan_firsts <- seq(as.Date(paste0(min(c(years, yearsU)), "-01-01")),
                  as.Date(paste0(max(c(years, yearsU)), "-01-01")),
                  by = "1 year")


p1 <- ggplot(uniprot_full.df, aes(x=Date, y=UniProtKB)) + 
  geom_area(color='black', fill='grey', size=1.5, show.legend = T) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = '', y= '') +
  scale_x_date(breaks = rev(rev(jan_firsts)[seq(1, length(jan_firsts), by = 2)]), date_labels = "%Y", limits = c(jan_firsts[1], jan_firsts[length(jan_firsts)]), expand = c(0,0)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())


p2 <- ggplot(genbank.df.l, aes(x=Date, y=entries, color=database, group = database)) + 
  # geom_line() +
  geom_area(aes(fill=database), alpha=0.3, size=1.5, position = 'identity') +
  scale_y_continuous(labels = label_comma()) +
  scale_x_date(breaks = rev(rev(jan_firsts)[seq(1, length(jan_firsts), by = 2)]), date_labels = "%Y", limits = c(jan_firsts[1], jan_firsts[length(jan_firsts)]), expand = c(0,0)) +
  theme_bw() + 
  labs(x = '', y= '') +
  scale_fill_manual(values = c('WGSSequences' = '#009bff', 'GenBankSequences' = 'black' )) +
  scale_color_manual(values = c('WGSSequences' = '#009bff', 'GenBankSequences' = 'black' )) +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

p2 / p1

