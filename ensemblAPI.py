#!/usr/bin/env python

# Python 3 compatible

import urllib.request
import urllib.parse
import urllib.error

import json
import time
import re
import sys


class EnsemblRestClient(object):
    def __init__(self, server='http://rest.ensembl.org/', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None, doseq=0):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urllib.parse.urlencode(params, doseq)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = urllib.request.Request(self.server + endpoint, headers=hdrs)
            response = urllib.request.urlopen(request)
            content = response.read().decode('utf-8')

            if content:
                data = json.loads(content)
            self.req_count += 1

        except urllib.error.HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry) + 1.0)
                    self.perform_rest_action(endpoint, hdrs)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def get_ensembl_id(self, species, symbol, only_stable_id=True):
        """
        Return an ensembl ID for the species gene symbol.

        Parameters:
        -----------
        species         -   A species name, e.g., human
        symbol          -   A gene symbol, e.g., BRAF
        only_stable_id  -   If true returns only the first gene in the list of genes
                            given by the sever.
                            else it will return the full JSON object.
        """

        genes = self.perform_rest_action(
            '/xrefs/symbol/{0}/{1}'.format(species, symbol),
            params={'object_type': 'gene'}
        )

        if genes:
            if only_stable_id:
                stable_id = genes[0]['id']
                return(stable_id)

            return(genes)

        return None

    def get_variants(self, species, symbol):
        stable_id = self.get_ensembl_id(species, symbol)

        if stable_id:
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'}
            )
            return variants

        return None

    def get_ensembl_info(self, identifier, expand=False):
        """
        Returns basic information on the given identifier.

        Parameters:
        -----------
        identifier  -   Ensembl identifier, e.g., ENSE00001939332.
        expand      -   Expands the search to include any connected features,
                                e.g., if the object is a gene, its transcripts,
                                translations and exons will be returned as well.
        """
        if expand:
            ext = "/lookup/id/{0}?expand=1".format(identifier)
        else:
            ext = "/lookup/id/{0}?".format(identifier)

        info = self.perform_rest_action(ext)

        if info:
            return(info)

        return(None)

    def get_overlap(self, species='human', region='7:140424943-140624564', features=None):
        """Return

        Parameters:
        -----------
        species     -   The organism whose genome is to be searched.
        region      -   The region within the genom to search.
        features    -   A list of feature to retrieve. Multiple values are accepted.
                        Enum(gene, transcript, cds, exon, repeat, simple, misc, variation,
                        somatic_variation, structural_variation, somatic_structural_variation,
                        constrained, regulatory, segmentation, motif, chipseq, array_probe).
        """
        features = features or ['gene']
        # Format the region into an appropriate format for query.
        region = self.format_region(region)

        info = self.perform_rest_action(
            '/overlap/region/{0}/{1}'.format(species, region),
            params={'feature': features},
            doseq=1
        )

        if info:
            return(info)

        return(None)

    @staticmethod
    def format_region(region):
        """Return the region in the correct format for submission to Ensembl.

            Parameters:
            -----------
            region  -   Either a coordinate string, chr7:123311-212131-,
                        or a list, ['7', '123311', '212131']

        """
        try:
            m = re.search(r'([XY]|[0-9]*):[0-9]*-[0-9]*', region)
            return(m.group())
        except TypeError:
            return(region[0] + ':' + str(region[1]) + '-' + str(region[2]))


if __name__ == '__main__':
    if len(sys.argv) == 3:
        species, symbol = sys.argv[1:]
    else:
        species, symbol = 'mouse', 'BRAF'

    # Optional server
    # client = EnsemblRestClient(server='http://grch37.rest.ensembl.org/')
    client = EnsemblRestClient()

    variants = client.get_variants(species, symbol)
    if variants:
        for v in variants:
            print('{seq_region_name}:{start}-{end}:{strand} ==> {id} ({consequence_type})'.format(**v))
    print()

    genes = client.get_ensembl_id(species, symbol)
    print(genes, end='\n\n')

    overlapping = client.get_overlap(species='human')
    print(overlapping, end='\n\n')
