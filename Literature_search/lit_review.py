import pandas as pd

pubmed = pd.read_csv("~/Dropbox/Literature_reviews/Burkitt_mutations/pubmed_export.csv")
wos = pd.read_csv("~/Dropbox/Literature_reviews/Burkitt_mutations/webofscience_export.csv")
scopus = pd.read_csv("~/Dropbox/Literature_reviews/Burkitt_mutations/scopus_export.csv")

print(pubmed.shape)
print(wos.shape)
print(scopus.shape)

all_results = pubmed.append(wos)
# all_results = all_results.append(scopus)
all_results["doi_lower"] = all_results["DOI"].str.lower()
all_results = all_results.drop_duplicates("doi_lower")

print(all_results.shape)
print(all_results.head())

all_results.to_csv("~/Dropbox/Literature_reviews/Burkitt_mutations/duplicates_removed.csv")
