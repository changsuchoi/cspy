# daily tns report
# unclassified AT
# classified object
# include FRB
# save csv file as tns_search.csv in /home/cschoi/Download
# https://www.wis-tns.org/content/tns-getting-started

import os
import astropy.io.ascii as ascii
import sys
from astropy.table import Table, Column
import astropy.coordinates as coord
from astropy.time import Time
import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.coordinates import ICRS
from astropy.coordinates import SkyCoord
import pandas as pd
from datetime import datetime

os.chdir('/data7/cschoi/sngal/recent-sne-check/rochester-list')
#datetime.today().strftime("%Y%m%d%H%M%S")    # YYYYmmddHHMMSS 형태의 시간 출력
#datetime.today().strftime("%Y/%m/%d %H:%M:%S")  # YYYY/mm/dd HH:MM:SS 형태의 시간 출력
today=datetime.today()
today=datetime.today().strftime("%Y%m%d %H:%M:%S")[:8]
print(today)

radius=10.0  # 30 arcmin = 0.5 deg

print ('Radius '+str(radius)+' arcmin')
# print ('Reading recentsnelist.txt file ...')
# colnames=['Ra','Dec','EarliestObs','Host','Type','Last','Max','Link','Discoverer']
# latestsnelist=pd.read_table('recentsnelist.txt')
# latestsnelist=pd.read_table('recentlist.txt')	#,names=colnames,data_start=1,guess='False')
# latestsnelist=ascii.read('recentsnelist.txt',delimiter='\t')	#,names=colnames,data_start=1,guess='False')
imsnglist=ascii.read('/data7/cschoi/IMSNG/target/alltarget.dat')

# tns server search for 7days
# three 500 num pages files
os.system("curl \'https://www.wis-tns.org/search?&page=0&discovered_period_value=7&discovered_period_units=days&unclassified_at=0&classified_sne=0&include_frb=1&name=&name_like=0&isTNS_AT=all&public=all&ra=&decl=&radius=&coords_unit=arcsec&reporting_groupid%5B%5D=null&groupid%5B%5D=null&classifier_groupid%5B%5D=null&objtype%5B%5D=null&at_type%5B%5D=null&date_start%5Bdate%5D=&date_end%5Bdate%5D=&discovery_mag_min=&discovery_mag_max=&internal_name=&discoverer=&classifier=&spectra_count=&redshift_min=&redshift_max=&hostname=&ext_catid=&ra_range_min=&ra_range_max=&decl_range_min=&decl_range_max=&discovery_instrument%5B%5D=null&classification_instrument%5B%5D=null&associated_groups%5B%5D=null&official_discovery=0&official_classification=0&at_rep_remarks=&class_rep_remarks=&frb_repeat=all&frb_repeater_of_objid=&frb_measured_redshift=0&frb_dm_range_min=&frb_dm_range_max=&frb_rm_range_min=&frb_rm_range_max=&frb_snr_range_min=&frb_snr_range_max=&frb_flux_range_min=&frb_flux_range_max=&num_page=500&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bremarks%5D=0&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0&display%5Brepeater_of_objid%5D=0&display%5Bdm%5D=0&display%5Bgalactic_max_dm%5D=0&display%5Bbarycentric_event_time%5D=0&display%5Bpublic_webpage%5D=0&format=csv\' > tns0.csv")
os.system("curl \'https://www.wis-tns.org/search?&page=1&discovered_period_value=7&discovered_period_units=days&unclassified_at=0&classified_sne=0&include_frb=1&name=&name_like=0&isTNS_AT=all&public=all&ra=&decl=&radius=&coords_unit=arcsec&reporting_groupid%5B%5D=null&groupid%5B%5D=null&classifier_groupid%5B%5D=null&objtype%5B%5D=null&at_type%5B%5D=null&date_start%5Bdate%5D=&date_end%5Bdate%5D=&discovery_mag_min=&discovery_mag_max=&internal_name=&discoverer=&classifier=&spectra_count=&redshift_min=&redshift_max=&hostname=&ext_catid=&ra_range_min=&ra_range_max=&decl_range_min=&decl_range_max=&discovery_instrument%5B%5D=null&classification_instrument%5B%5D=null&associated_groups%5B%5D=null&official_discovery=0&official_classification=0&at_rep_remarks=&class_rep_remarks=&frb_repeat=all&frb_repeater_of_objid=&frb_measured_redshift=0&frb_dm_range_min=&frb_dm_range_max=&frb_rm_range_min=&frb_rm_range_max=&frb_snr_range_min=&frb_snr_range_max=&frb_flux_range_min=&frb_flux_range_max=&num_page=500&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bremarks%5D=0&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0&display%5Brepeater_of_objid%5D=0&display%5Bdm%5D=0&display%5Bgalactic_max_dm%5D=0&display%5Bbarycentric_event_time%5D=0&display%5Bpublic_webpage%5D=0&format=csv\' > tns1.csv")
os.system("curl \'https://www.wis-tns.org/search?&page=2&discovered_period_value=7&discovered_period_units=days&unclassified_at=0&classified_sne=0&include_frb=1&name=&name_like=0&isTNS_AT=all&public=all&ra=&decl=&radius=&coords_unit=arcsec&reporting_groupid%5B%5D=null&groupid%5B%5D=null&classifier_groupid%5B%5D=null&objtype%5B%5D=null&at_type%5B%5D=null&date_start%5Bdate%5D=&date_end%5Bdate%5D=&discovery_mag_min=&discovery_mag_max=&internal_name=&discoverer=&classifier=&spectra_count=&redshift_min=&redshift_max=&hostname=&ext_catid=&ra_range_min=&ra_range_max=&decl_range_min=&decl_range_max=&discovery_instrument%5B%5D=null&classification_instrument%5B%5D=null&associated_groups%5B%5D=null&official_discovery=0&official_classification=0&at_rep_remarks=&class_rep_remarks=&frb_repeat=all&frb_repeater_of_objid=&frb_measured_redshift=0&frb_dm_range_min=&frb_dm_range_max=&frb_rm_range_min=&frb_rm_range_max=&frb_snr_range_min=&frb_snr_range_max=&frb_flux_range_min=&frb_flux_range_max=&num_page=500&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bremarks%5D=0&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0&display%5Brepeater_of_objid%5D=0&display%5Bdm%5D=0&display%5Bgalactic_max_dm%5D=0&display%5Bbarycentric_event_time%5D=0&display%5Bpublic_webpage%5D=0&format=csv\' > tns2.csv")

# one day search
#os.system("curl \'https://www.wis-tns.org/&page=X&search?&discovered_period_value=7&discovered_period_units=days&unclassified_at=0&classified_sne=0&include_frb=1&name=&name_like=0&isTNS_AT=all&public=all&ra=&decl=&radius=&coords_unit=arcsec&reporting_groupid%5B%5D=null&groupid%5B%5D=null&classifier_groupid%5B%5D=null&objtype%5B%5D=null&at_type%5B%5D=null&date_start%5Bdate%5D=&date_end%5Bdate%5D=&discovery_mag_min=&discovery_mag_max=&internal_name=&discoverer=&classifier=&spectra_count=&redshift_min=&redshift_max=&hostname=&ext_catid=&ra_range_min=&ra_range_max=&decl_range_min=&decl_range_max=&discovery_instrument%5B%5D=null&classification_instrument%5B%5D=null&associated_groups%5B%5D=null&official_discovery=0&official_classification=0&at_rep_remarks=&class_rep_remarks=&frb_repeat=all&frb_repeater_of_objid=&frb_measured_redshift=0&frb_dm_range_min=&frb_dm_range_max=&frb_rm_range_min=&frb_rm_range_max=&frb_snr_range_min=&frb_snr_range_max=&frb_flux_range_min=&frb_flux_range_max=&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bremarks%5D=0&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0&display%5Brepeater_of_objid%5D=0&display%5Bdm%5D=0&display%5Bgalactic_max_dm%5D=0&display%5Bbarycentric_event_time%5D=0&display%5Bpublic_webpage%5D=0&format=csv\'  > tnstest.csv")

tns0=ascii.read('tns0.csv')
tns1=ascii.read('tns1.csv')
tns2=ascii.read('tns2.csv')

# vstack
import

































'''
https://www.wis-tns.org/search?
&page=0
&discovered_period_value=1
&discovered_period_units=days
&unclassified_at=0
&classified_sne=0
&include_frb=1
&name=
&name_like=0
&isTNS_AT=all
&public=all
&ra=
&decl=
&radius=
&coords_unit=arcsec
&reporting_groupid%5B%5D=null
&groupid%5B%5D=null
&classifier_groupid%5B%5D=null
&objtype%5B%5D=null
&at_type%5B%5D=null
&date_start%5Bdate%5D=
&date_end%5Bdate%5D=
&discovery_mag_min=
&discovery_mag_max=
&internal_name=
&discoverer=
&classifier=
&spectra_count=
&redshift_min=
&redshift_max=
&hostname=
&ext_catid=
&ra_range_min=
&ra_range_max=
&decl_range_min=
&decl_range_max=
&discovery_instrument%5B%5D=null
&classification_instrument%5B%5D=null
&associated_groups%5B%5D=null
&official_discovery=0
&official_classification=0
&at_rep_remarks=
&class_rep_remarks=
&frb_repeat=all
&frb_repeater_of_objid=
&frb_measured_redshift=0
&frb_dm_range_min=
&frb_dm_range_max=
&frb_rm_range_min=
&frb_rm_range_max=
&frb_snr_range_min=
&frb_snr_range_max=
&frb_flux_range_min=
&frb_flux_range_max=
&num_page=500
&display%5Bredshift%5D=1
&display%5Bhostname%5D=1
&display%5Bhost_redshift%5D=1
&display%5Bsource_group_name%5D=1
&display%5Bclassifying_source_group_name%5D=1
&display%5Bdiscovering_instrument_name%5D=0
&display%5Bclassifing_instrument_name%5D=0
&display%5Bprograms_name%5D=0
&display%5Binternal_name%5D=1
&display%5BisTNS_AT%5D=0
&display%5Bpublic%5D=1
&display%5Bend_pop_period%5D=0
&display%5Bspectra_count%5D=1
&display%5Bdiscoverymag%5D=1
&display%5Bdiscmagfilter%5D=1
&display%5Bdiscoverydate%5D=1
&display%5Bdiscoverer%5D=1
&display%5Bremarks%5D=0
&display%5Bsources%5D=0
&display%5Bbibcode%5D=0
&display%5Bext_catalogs%5D=0
&display%5Brepeater_of_objid%5D=0
&display%5Bdm%5D=0
&display%5Bgalactic_max_dm%5D=0
&display%5Bbarycentric_event_time%5D=0
&display%5Bpublic_webpage%5D=0
&page=X
&format=csv



Quick query links
Unclassified ATs:

https://www.wis-tns.org/search?&isTNS_AT=yes&unclassified_at=1&num_page=500

Unclassified ATs within a discovery magnitude range (e.g. 15-17):

https://www.wis-tns.org/search?&isTNS_AT=yes&unclassified_at=1&discovery...

Classified SNe:

https://www.wis-tns.org/search?&isTNS_AT=yes&classified_sne=1&num_page=500



To download as CSV, add "&format=csv", e.g. for classified SNe:

https://www.wis-tns.org/search?&isTNS_AT=yes&classified_sne=1&num_page=5...

Additional notes:

- The maximum allowed number of results per page is 500 (num_page=500).

- If downloading as CSV/TSV and the retrieved results surpass the maximum, you
may implement paging by specifying &page=X (0 being the 1st page), e.g.:

https://www.wis-tns.org/search?&page=0&isTNS_AT=yes&classified_sne=1&num...

https://www.wis-tns.org/search?&page=1&isTNS_AT=yes&classified_sne=1&num...

[ Top ]




Daily CSV Staging

Every day after UT midnight, two CSV files are created and are accessible for
download under: https://www.wis-tns.org/system/files/tns_public_objects/

1. tns_public_objects.csv.zip - holds the entire catalog of TNS public objects
(AT/SN/FRB/... ~70,000 currently). This file is overwritten daily. The date and
time when the list of objects was created is specified in the first line; e.g.
"2021-03-15 00:00:00"

2. tns_public_objects_YYYYMMDD.csv.zip - holds only those entries (objects) that
were either added or modified during the specified day. So, e.g. during Mar 15,
2021 it is possible to download this latest CSV for the previous day:
tns_public_objects_20210314.csv.zip The first line in the CSV will contain the
exact duration covering the entries in the file; e.g. for the above example:
"2021-03-14 00:00:00 - 23:59:59"

The separate daily files remain in place for 1 month backwards.

Staging the CSV files serves to fulfil requests by TNS users, as well as encourage
performing time-consuming operations locally by users, reducing the load on the
TNS servers.

For example, if you need to cross-match entire catalogs or long object lists, we
request that this would be done locally, against the csv (or a locally managed
DB), rather than by executing multiple cone-searches via the Search API. Calling
the APIs for a limited number of objects is clearly fine, but we ask that our
users apply appropriate caution and sensibility when using the TNS resources,
that serve a broad community.

The csv's contain the following columns:

"objid","name_prefix","name","ra","declination","redshift","typeid","type",
"reporting_groupid","reporting_group","source_groupid","source_group",
"discoverydate","discoverymag","discmagfilter","filter",
"reporters","time_received","internal_names","creationdate","lastmodified"
cURL usage examples for downloading (note that the api_key is required):

curl -X POST -d 'api_key=YOUR-API-KEY' https://www.wis-tns.org/system/files/tns_public_objects/tns_public_objec... > tns_public_objects.csv.zip
curl -X POST -d 'api_key=YOUR-API-KEY' https://www.wis-tns.org/system/files/tns_public_objects/tns_public_objec... > tns_public_objects_20210314.csv.zip
Or if logged in, you can just go to the URL for downloading the CSV:

https://www.wis-tns.org/system/files/tns_public_objects/tns_public_objec...
https://www.wis-tns.org/system/files/tns_public_objects/tns_public_objec... (specify the required date)
'''
