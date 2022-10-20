import os
import time
import sys
import argparse
import calendar
import requests
import numpy as np
from tqdm import tqdm
import json
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u

###

def get_token(username = None, password = None, token = None):
    
    if (username is None or password is None) and token is None:
        raise NameError('Either username AND password, or token must be given')
    elif (username is None or password is None) and token is not None:
        token = token
    elif (username is not None and password is not None) and token is None:
        try:
            response = requests.post('https://observe.lco.global/api/api-token-auth/',
                                     data = {'username': username,
                                             'password': password}).json()
            token = response['token']
        except:
            print('Hmm. That didnt work. Username/password probably not accepted. Proceeding as anonymous user.')
    elif username is not None and password is not None and token is not None:
        try:
            response = requests.post('https://observe.lco.global/api/api-token-auth/',
                                     data = {'username': username,
                                             'password': password}).json()
            token_retrieved = response['token']
            if token != token_retrieved:
                print('Retrieved token from username/password does not match given token. Using the token retrieved using the given username/password.')
        except:
            print('Username, password, AND token all givenm, but username/password not accepted. Using the given token only.')

    return token
        

def get_url(proposal_id,
            start,
            end,
            reduction_level,
            primary_optical_element,
            telescope_id,
            target_name,
            covers,
            limit = 1000, offset = 0):
    
    url = 'https://archive-api.lco.global/frames/?'\
        f'proposal_id={proposal_id}&'\
        f'start={start}&'\
        f'end={end}&'\
        f'covers={covers}&'\
        'public=true&'\
        f'limit={limit}&'\
        f'offset={offset}&'\
        f'target_name={target_name}&'\
        f'reduction_level={reduction_level}&'\
        f'primary_optical_element={primary_optical_element}&'\
        f'telescope_id={telescope_id}'
    
    return url

def get_ra_dec(target_name):
    
    result_table = Simbad.query_object(target_name)
    # print(':'.join(result_table['RA'][0]), ':'.join(result_table['DEC'][0]))
    coord = SkyCoord(ra = result_table['RA'][0], dec = result_table['DEC'][0],
                     frame = 'icrs', unit = (u.hourangle, u.deg))
    ra, dec = coord.ra.deg, coord.dec.deg

    return ra, dec

def get_frames(url, username, password, token):
    
    print(f'Fetching frames for the url: {url}')
    # print(f'{username=} {password=} {token=}')
    
    token = get_token(username, password, token)
    if token is None:
        if (username is not None or password is not None):
            print('Token could not be retrieved. Proceeding as anonymous user.')
        response = requests.get(url).text
    else:
        response = requests.get(url, headers = {'Authorization' : f'token {token}'}).text
    
    response = json.loads(response)
    frames = response['results']
    if response['next']:
        more = True
        response_new = response
        while more:
            response_new = requests.get(response_new['next'], headers = {'Authorization' : f'token {token}'}).json()
            frames += response_new['results']
            if response_new['next']:
                more = True
            else:
                more = False
    
    return frames

def download_frames(frames, n_chunk = 100):
    
    if len(frames) > n_chunk:
        return

def parse_args():
    
    """Parse command-line inputs"""

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-start', '--start', default='',
                        help='Start date for search [YYYY-MM-DD]')
    parser.add_argument('-end', '--end', default='',
                        help='End date for search. [YYYY-MM-DD]')
    
    parser.add_argument('-id', '--proposal_id', default='',
                        help='List of proposals to search for')
    parser.add_argument('-targ', '--target_name', default='',
                        help='Object name to search for')
    parser.add_argument('-red_lev', '--reduction_level', default='91',
                        help='Reduction level to filter by. 91: BANZAI')
    parser.add_argument('-filters', '--filters', default='',
                        help='Filters to get data for. Default is to collect images in all filters.', type = str)
    parser.add_argument('-telescope_ids', '--telescope_ids', default = '',
                        help = 'Telescope ID to restrict search to. \'0m4a\', \'1m0a\', and \'2m0a\' are the 0.4-, 1-, and 2-meter telescopes, respectively.')
    
    parser.add_argument('-lim', '--limit', default='1000',
                        help='Number of frames to limit the search to')
    
    parser.add_argument('-outdir', '--outdir', default='./',
                        help='Directory where the data will be downloaded into. Default is the current directory')
    parser.add_argument('--overwrite', action = 'store_true', default = False,
                        help='If true, overwrite files even if they already exist in \'outfir\'') 
    parser.add_argument('-token', '--token', default = None,
                        help = 'Authorization token to use. Default is that of John Wendeborn')
    parser.add_argument('-username', '--username', default = None,
                        help = 'Username for token authorization')
    parser.add_argument('-password', '--password', default = None,
                        help = 'Password for authorization token')
    parser.add_argument('-verbose', '--verbose', default = False, action = 'store_true',
                        help = '')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    
    args = parse_args()
    
    proposal_id = args.proposal_id
    start = args.start
    end = args.end
    
    target_name = args.target_name
    reduction_level = args.reduction_level
    filters = args.filters
    if filters in ['all', '', '[all]', '[\'all\']']:
        filters = ['']
    else:
        if ',' in filters:
            filters = filters.split(',')
        else:
            filters = [filters]
    telescope_ids = args.telescope_ids
    if telescope_ids in ['all', '']:
        telescope_ids = ['']
    else:
        if ',' in telescope_ids:
            telescope_ids = telescope_ids.split(',')
        else:
            telescope_ids = [telescope_ids]
    
    limit = args.limit
    
    username = args.username
    password = args.password
    token = args.token
    
    outdir = args.outdir
    overwrite = args.overwrite
    verbose = args.verbose
    
    covers = ''
    if target_name != '':
        try:
            ra, dec = get_ra_dec(target_name = target_name)
            covers = f'POINT%28{ra}%20{dec}%29'
            target_name = ''
        except:
            pass
    
    all_frames = []
    for filt in filters:
        for telescope_id in telescope_ids:
            request_url = get_url(proposal_id = proposal_id, start = start, end = end, 
                                  reduction_level = reduction_level, primary_optical_element = filt,
                                  target_name = target_name, telescope_id = telescope_id,
                                  limit = limit, covers = covers)
            
            frames = get_frames(request_url, username = username, password = password, token = token)
            all_frames += frames
    
    existing_files, new_files, failed_files = 0, 0, 0
    for frame in tqdm(all_frames, desc = f'Found {len(all_frames)} frames. Downloading to {outdir}'):
        url, filename = frame['url'], frame['filename']
        if verbose:
            print(filename)
        outfile = os.path.join(outdir, filename)
        if os.path.exists(outfile):
            if overwrite:
                if verbose:
                    print(f'{filename} already exists. Overwriting.')
            elif os.path.getsize(outfile) < 1:
                if verbose:
                    print('File exists, but it is 0 bytes. Overwriting.')
            else:
                if verbose:
                    print(f'{filename} already exists. Will not overwrite.')
                existing_files += 1
                continue
        with open(outfile, 'wb+') as f:
            # try:
            f.write(requests.get(url).content)
            new_files += 1
            # except:
            #     failed_files += 1
                
    print(f'Total Files Found: {len(frames)}')
    print(f'    {new_files} New')
    print(f'    {existing_files} Already Exist')
    print(f'    {failed_files} Failed to Download')


