import requests
headers = {'Authorization': 'sturkarslan@gmail.com 20444f8b27f4c60aaba2'}
payload = {'cancer_type': 'NS', 'title': 'Test', 'reference': 'hg38'}
r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                headers=headers,
                files={
                        'mutations': open('P76156_6_CGI_input.txt', 'rb'),
                        },
                data=payload)
r.json()