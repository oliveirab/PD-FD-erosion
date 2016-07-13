FindIUCN <- function(input, count = FALSE) {
  
  input <- b.getnames(input)
  
  # Matrices to save the result
  ln <- length(input)
  matriz1 <- matrix(nrow = ln, ncol = 1)
  IUCNname <- matriz1
  status <-  matriz1
  criterio <- matriz1
  populacao <- matriz1
  familia <- matriz1
  autor <- matriz1  
  pais <- matriz1
  
  # With count window
  if (count) {
    
    # Do not set a new device in rstudio to avoid warnings()
    if (!"tools:rstudio" %in% search()) {
      dev.new(width = 2, height = 2, pointsize = 12)
      par(mar = c(0, 0, 0, 0))
    }
    
    for (i in 1:ln){
      plot.new()
      text(0.5, 0.5, paste(paste("Total:", ln, "\n",
                                 "Species to go: ",
                                 (ln - i))))
      
      loopresu        <- b.LoopIUCN(input[i])
      IUCNname[i, ]    <- loopresu$IUCNname
      status[i, ]    <- loopresu$Status
      criterio[i, ]  <- loopresu$Criteria
      populacao[i, ] <- loopresu$Population
      familia[i, ]   <- loopresu$Family
      autor[i, ]     <- loopresu$Author
      pais[i, ]       <- loopresu$Country
    }
    dev.off()
  }
  
  if(!count){
    for (i in 1:ln){          
      
      loopresu        <- b.LoopIUCN(input[i])
      IUCNname[i, ]    <- loopresu$IUCNname
      status[i, ]    <- loopresu$Status
      criterio[i, ]  <- loopresu$Criteria
      populacao[i, ] <- loopresu$Population
      familia[i, ]   <- loopresu$Family
      autor[i, ]     <- loopresu$Author
      pais[i, ]       <- loopresu$Country
    }
  }
  
  # Final table
  resu <- cbind(input, IUCNname, familia, status, criterio, populacao, autor, pais)
  colnames(resu) <- c("Species", "IUCNname", "Family", "Status", "Criteria",
                      "Population", "Description_Year", "Country")
  
  # Making a data frame and guarantee that it is numeric the Year
  resu <- as.data.frame(resu)
  resu[, 6] <- as.numeric(levels(resu[, 6]))[resu[, 6]]
  return(resu)
}





# Inside the loop function -----------------------------

b.LoopIUCN <- function(input) {
  library(XML)
  
  # Get species code
  speciescode <- b.getcode(input)
  # Find the web
  h <- try(htmlParse(paste("http://api.iucnredlist.org/details/",
                           speciescode, "/0", sep = "")),
           silent=TRUE)
  # If not find, species set to NE
  if (class(h)[1] == "try-error") {
    name <- "NA"
    statusi <- "NE"
    criterioi <- ""
    populacaoi <- "Unknown"
    familiai <- ""
    autori <- ""
    paisi <- ""
  } else {
    # Look for the informations
    name <- try(xpathSApply(h, '//h1[@id="scientific_name"]', 
                                xmlValue),
                silent = TRUE)
    statusi <- try(xpathSApply(h, '//div[@id="red_list_category_code"]', 
                               xmlValue),
                   silent = TRUE)
    criterioi <- try(xpathSApply(h, '//div[@id="red_list_criteria"]',
                                 xmlValue),
                     silent = TRUE)
    pop <- try(xpathSApply(h, '//div[@id="population_trend"]',
                           xmlValue),
               silent=TRUE)
    
    
    familiai <- try(xpathSApply(h, '//div[@id="family"]',
                                xmlValue),
                    silent = TRUE)
    
    autori <- try(xpathSApply(h, '//div[@id="species_authority"]',
                              xmlValue),
                  silent = TRUE)
    
    # Error control
    if (class(name)[1] == "try-error") {
      name <-  "NA"
    } else {
      name <- name
    }
    
    if (class(statusi)[1] == "try-error") {
      statusi <-  "NE"
    } else {
      statusi <- statusi
    }
    
    if (class(criterioi)[1] == "try-error") {
      criterioi <-  ""
    } else {
      criterioi <- criterioi
    }
    
    # Sometimes the pop is an empty list.
    if (class(pop)[1] == "try-error" | class(pop) == "list") {
      populacaoi <-  "Unknown"
    } else {
      populacaoi <- pop
    }
    
    if (class(familiai)[1] == "try-error") {
      familiai <-  ""
    } else {
      familiai <- familiai
    }
    
    if (class(autori)[1] == "try-error") {
      autori <-  ""
    } else {
      autori <- gsub("\\D", "", autori)
      autori <- as.numeric(substr(autori, 1, 4))
      if (is.na(autori)) {
        autori <- ""
      }
    }
    
    # Countries
    distr1 <- try(xpathSApply(h, '//ul[@class="countries"]', 
                              xmlValue),
                  silent = TRUE)
    # Error control
    if (class(distr1)[1] == "try-error") {
      paisi <- ""
    } else {
      # Remove countries with comma
      distr2 <- try(unlist(strsplit(distr1, "\n")), silent = TRUE)
      distr2 <- gsub(',', '', distr2)
      paisi <- paste(distr2, collapse = ", ")
    }
  }
  # Return the information
  return(list("IUCNname" = name,
              "Status" = statusi,
              "Criteria" = criterioi,
              "Population" = populacaoi,
              "Family" = familiai,
              "Author" = autori,
              "Country" = paisi))
}

# Inside the loop function -----------------------------

b.getcode <- function(input) {
  
  binomialerror <- length(unlist(strsplit(input, " "))) == 2
  input <- gsub(as.matrix(input), pattern = " ", replacement = "-")
  
  h <- try(htmlParse(paste("http://api.iucnredlist.org/go/",
                           input, sep = "")),
           silent = TRUE)
  if (class(h)[1] != "try-error" & binomialerror) {
    b <- try(xpathSApply(h, '//div', xmlValue), silent = TRUE)[1]
    c <- as.numeric(gsub("\\D", "", b))
    
    
    # Subsecies control
    http <- "http://www.iucnredlist.org/details/summary/"
    h1 <- htmlParse(paste(http, c, "/0", sep = ""))
    links <- xpathSApply(h1, "//a/@href")
    links <- strsplit(links, "\n")
    parents <- xpathSApply(h1, "//a")
    
    # function to transform xml in character list
    tocharacter <- function(x) {
      do.call(paste, as.list(capture.output(x)))
    }
    parents2 <- sapply(parents, tocharacter)
    #menos <- length(parents2) - length(links)
    
    posParent <- grep("_parent", parents2)
    if (length(posParent) == 1) {
      #cpar <- gsub("\\D", "", (links[posParent - menos]))
      cpar <- gsub("\\D", "", (parents2[posParent]))
      c <- as.numeric(substr(cpar, 1, nchar(cpar) - 1))
    }
    ################################################
    return(c)
  } else {
    return(NULL)
  }
}

# Inside the loop function -----------------------------
b.getnames <- function(input) {
  
  # Get species from a PAM
  if (class(input) == "PresenceAbsence") {
    input <- input$S
  }
  
  # Accept species separeted by underline or space
  input <- gsub(as.matrix(input), pattern = "_", 
                replacement = " ")
  
  # Remove double or more spaces
  input <- gsub("\\s{2, }", " ", input)
  
  # Remove space from the beggining and end
  trim <- function(x) { 
    return(gsub("^\\s+|\\s+$", "", x))
  }
  
  input <- as.vector(trim(input))
  
  # Species with wrong names
  count2 <- function(x) {
    return(length(x) != 2)
  }
  
  binomialerror <- sapply((strsplit(input, " ")), count2)
  sps <- which(binomialerror)
  sps_name <- paste("\t", input[sps], "\n")
  
  # Error in species name control
  if (length(sps) > 0) {
    warning(paste("The following species do not follow a binomial nomeclature:\n",
                  paste(sps_name, collapse = "")))
  }
  return(input)
}