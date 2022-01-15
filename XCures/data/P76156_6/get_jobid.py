import requests
headers = {'Authorization': 'sturkarslan@gmail.com 20444f8b27f4c60aaba2'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
r.json()
print(r.json)