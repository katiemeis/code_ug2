#!/afs/crc.nd.edu/x86_64_linux/p/python/3.5.2/build/bin/python3
import sys
import os
import re
from lxml import html
from lxml import etree
import xml.etree.ElementTree as ET
import requests
import subprocess
#======================================================================
# Functions
#======================================================================
def scrap_source(R_url):
    html_data = requests.get(R_url)
    root = html.fromstring(html_data.content)
    ob = root.findall(".//td")
    package_name = re.findall("<h2>(.*?):.*?</h2>",html_data.text,re.IGNORECASE|re.MULTILINE|re.DOTALL)[0]
    source_html = ''
    for n in ob:
        a = html.tostring(n,encoding='utf8',method = 'text').decode('utf8')
        if re.search('Package',a):
            source_html = n.getnext()
    if len(source_html):
        package_file = html.tostring(source_html,encoding='utf8',method='text').decode('utf8').strip()

        return([package_file,R_url + '/' + source_html.xpath('./a/@href')[0],package_name])
    else:
        return(False)

def scrap_imports(R_url):
    html_data = requests.get(R_url)
    root = html.fromstring(html_data.content)
    ob = root.findall(".//td")
    imports = []
    for n in ob:
        a = html.tostring(n,encoding='utf8',method = 'text').decode('utf8')
        if re.search('Imports:',a):
            imports_html = n.getnext()
            imports_tmp = imports_html.xpath('./a/@href')
            for i in imports_tmp:                
                imports.append(R_url + '/' + i.replace('/index.html',''))
        if re.search('Depends:',a):
            imports_html = n.getnext()
            imports_tmp = imports_html.xpath('./a/@href')
            for i in imports_tmp:                
                imports.append(R_url + '/' + i.replace('/index.html',''))
    return(imports)

def check_package(p_name):
    # first check if package is installed
    str_run = "Rscript -e \'library(%s)\'" % (p_name)
    lib_status = subprocess.run(str_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
    lib_error = lib_status.stderr.decode()
    if re.search('error in library',lib_error,re.IGNORECASE):
        print('library %s not installed... proceeding to install with dependencies' % (p_name))
        return False
    else:
        print('library %s already installed and it seems to be working' % (p_name))
        return True

def install_package(r_url,l_dir,p_name):
    is_installed = check_package(p_name)
    if is_installed:
        return
    print(r_url)
    url_source = scrap_source(r_url)
    if not url_source:
        print('Invalid R package: %s\n' % (p_name),end='')
        quit()
    imports_list = scrap_imports(r_url)
    print(imports_list)
    for i in imports_list:
        print(i)
        source_import = scrap_source(i)
        install_package(i,l_dir,source_import[-1])

    # Finally install the package
    os.system("curl %s > %s" % (url_source[1], url_source[0]))
    os.system("R CMD INSTALL -l %s %s" % (l_dir,url_source[0]))
    os.system("rm %s" % (url_source[0]))

def main():    
    if len(sys.argv) < 3:
        print('''./install_R_packages.py [libdir] [package name]
        please spacify the filename and library directory''')
        quit()
    libdir = sys.argv[1]
    package_name = sys.argv[2]
    url_base = "https://cran.r-project.org/web/packages/";
    url_package = url_base + package_name
    install_package(url_package,libdir,package_name)
    os.system("fs setacl %s nd_campus rlidw" % (libdir))
    os.chdir(libdir)
    os.system("find . -type d -exec fs sa {} nd_campus read \;")


    
#======================================================================
# Generate submission files
#======================================================================
if __name__ == '__main__':
    main()
