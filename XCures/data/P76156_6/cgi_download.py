                    
import requests
headers = {'Authorization': 'sturkarslan@gmail.com 20444f8b27f4c60aaba2'}
payload={'action':'download'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/144a37d94b1c5064fb2f', headers=headers, params=payload)
with open('file.zip', 'wb') as fd:
    fd.write(r._content)