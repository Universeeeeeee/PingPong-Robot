"""Independent lower-level checks for raw morphology descriptors, not semantic labels."""
import hashlib
import json
import numpy as np
import review

OUT=review.OUT


def read_csv(name):
    import csv
    rows=list(csv.DictReader((OUT/name).open()))
    for r in rows:
        for k,v in r.items():
            if v in ['True','False']:r[k]=v=='True'
            elif not v:r[k]=np.nan
            else:
                try:r[k]=float(v)
                except ValueError:pass
    return rows


def main():
    rows=read_csv('morphology.csv');cards=read_csv('card_windows.csv');meta=json.loads((OUT/'review_metadata.json').read_text())
    assert len(rows)==630==len({(r['segment'],r['feature'],r['duration_ms'],r['frame']) for r in rows})
    assert len(cards)==90
    assert review.digest(review.B0)==meta['previous_B0_sha256']
    assert hashlib.sha256((OUT/'review.py').read_bytes()).hexdigest()==meta['review_script_sha256']
    assert hashlib.sha256((OUT/'protocol.md').read_bytes()).hexdigest()==meta['protocol_sha256']
    old={(r['segment'],r['feature'],r['duration_ms'],r['frame']):r for r in review.read_rows()}
    events,_=review.b0.raw.load_events();times=events['t'].astype(np.int64);t0=int(times[0]);cache={};checks=[]
    for r in rows:
        key=(r['segment'],r['feature'],r['duration_ms'],r['frame']);prior=old[key]
        assert r['position_valid']==prior['position_valid'] and r['status']==prior['status']
        if not r['position_valid']:
            assert not np.isfinite(r['raw_patch_events'])
            continue
        segment,duration,frame=r['segment'],int(r['duration_ms']),int(r['frame']);window_key=(segment,duration,frame)
        if window_key not in cache:
            lo=t0+round(review.b0.raw.SEGMENTS[segment][0]*1e6)+frame*duration*1000;hi=lo+duration*1000
            cache[window_key]=events[(times>=lo)&(times<hi)]
        e=cache[window_key];cx,cy=r['anchor_x'],r['anchor_y']
        local=e[(e['x'].astype(float)>=cx-13)&(e['x'].astype(float)<=cx+13)&(e['y'].astype(float)>=cy-13)&(e['y'].astype(float)<=cy+13)]
        points=set(zip(local['x'].tolist(),local['y'].tolist()))
        boundary=sum(abs(x-cx)>=12 or abs(y-cy)>=12 for x,y in points)
        # Explicit set flood-fill is independent of scipy.ndimage.label and raster indexing.
        remaining=points.copy();sizes=[]
        while remaining:
            stack=[remaining.pop()];size=0
            while stack:
                x,y=stack.pop();size+=1
                near={(x+dx,y+dy) for dx in [-1,0,1] for dy in [-1,0,1]}&remaining
                remaining-=near;stack.extend(near)
            sizes.append(size)
        assert len(local)==r['raw_patch_events']==prior['local_events']
        assert len(points)==r['distinct_pixels']==prior['local_distinct_pixels']
        assert abs(boundary/len(points)-r['boundary_2px_pixel_fraction'])<1e-12
        assert sum(size>=3 for size in sizes)==r['components_ge3_pixels']
        assert abs(max(sizes)/len(points)-r['largest_component_pixel_fraction'])<1e-12
        checks.append(dict(segment=segment,feature=r['feature'],duration_ms=duration,frame=frame,
                           event_payload_and_unique_pixels_match=True,independent_boundary_and_floodfill_match=True))
    for c in cards:
        key=(c['segment'],c['feature'],5.,c['review_5ms_frame']);r=old[key]
        if c['crop_source']=='valid_5ms_anchor':
            assert r['position_valid'] and c['crop_center_x']==r['anchor_x'] and c['crop_center_y']==r['anchor_y']
        else:
            assert not r['position_valid'] and [c['crop_center_x'],c['crop_center_y']]==review.b0.SEEDS[c['segment']][c['feature']]
        assert abs(c['midpoint_difference_2_minus_5_ms'])<=.50000001
    validation=dict(morphology_rows=630,card_panels=90,unique_keys=True,old_B0_files_unchanged=True,
                    valid_morphology_rows_independently_checked=len(checks),invalid_positions_not_filled=True,
                    descriptor_validation='Absolute timestamp mask; xy tuple sets; explicit 8-neighbor flood fill; border coordinates',
                    semantics_validated=False,semantics_note='Morphology descriptor implementation checks do not validate logo/Marker identity',
                    script_sha256=hashlib.sha256((OUT/'validate_review.py').read_bytes()).hexdigest(),
                    checks=checks)
    (OUT/'review_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps({k:v for k,v in validation.items() if k!='checks'},indent=2))


if __name__=='__main__':main()
