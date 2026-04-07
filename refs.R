library(RefManageR)

bib <- ReadBib("/Users/kellyclintdale/Documents/Research_Admin/kelly lab management/kelly-cv/kellyrefs.bib")

# Journal articles
articles <- bib[bib$bibtype == "Article"]

# Print nicely
print(articles, .opts = list(style = "text"))
