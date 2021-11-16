from bs4 import BeautifulSoup
import urllib.request
import json
import csv

user_agent = 'Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.0.7) Gecko/2009021910 Firefox/3.0.7'

url = "https://zoonomiaproject.org/the-mammal-tree-list-view/"
headers={'User-Agent':user_agent,} 

request=urllib.request.Request(url,None,headers) #The assembled request
response = urllib.request.urlopen(request)
data = response.read()
mammalian_tree = BeautifulSoup(data, "html.parser")

tree = mammalian_tree.find_all("div", {"class": "tree"})


with open('mammalian_web_scraped.csv', 'w') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
        filewriter.writerow(["SPECIES","ORDER","FAMILY","COMMON NAME"])

for t in tree:
    _species = str(t.find(class_ = "species").text.strip())
    _order = str(t.find(class_ = "order").text.strip())
    _family = str(t.find(class_ = "family").text.strip())
    _name = str(t.find(class_ = "common-name").text.strip())
    with open('mammalian.csv', 'a+', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
        filewriter.writerow([_species,_order,_family,_name])