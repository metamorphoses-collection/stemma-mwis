"""Transkribus Processing API helper. Credentials from the environment file named by TRANSKRIBUS_ENV (never printed)."""
import os, json, urllib.request, urllib.parse, time
def env():
    d = {}
    for l in open(os.environ['TRANSKRIBUS_ENV']):
        if '=' in l: k, v = l.strip().split('=', 1); d[k] = v.strip().strip('"').strip("'")
    return d
def token():
    e = env()
    data = urllib.parse.urlencode(dict(grant_type='password', client_id='processing-api-client',
                                       username=e['TRANSKRIBUS_USER'], password=e['TRANSKRIBUS_PASS'])).encode()
    r = urllib.request.urlopen(urllib.request.Request(
        'https://account.readcoop.eu/auth/realms/readcoop/protocol/openid-connect/token', data=data), timeout=60)
    return json.load(r)['access_token']
def api(tok, method, path, body=None):
    req = urllib.request.Request('https://transkribus.eu/processing/v1/' + path, method=method,
                                 data=json.dumps(body).encode() if body is not None else None,
                                 headers={'Authorization': 'Bearer ' + tok, 'Content-Type': 'application/json'})
    return json.load(urllib.request.urlopen(req, timeout=120))
if __name__ == '__main__':
    t = token(); print('token ok, length', len(t))
