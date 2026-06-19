import requests
import xml.etree.ElementTree as ET
import csv

import MySQLdb as mydb

import time
import datetime

from lxml import etree
import re

import pandas as pd
import numpy as np

from urllib.request import urlopen

import sys 

import http.client



"""
uid = "P0DTC2"
# UniProt XMLデータのURL
url = "https://rest.uniprot.org/uniprotkb/{}.xml".format(uid)
response = requests.get(url)
response.raise_for_status()

# XMLをパース
root = ET.fromstring(response.content)

# 名前空間を取得（UniProt XMLは名前空間付き）
ns = {'u': 'http://uniprot.org/uniprot'}

# 対象とするtype
target_types = {"domain", "region of interest", "topological domain", "site", "transmembrane region"}

# 結果格納用リスト
features = []

for feature in root.findall(".//u:feature", ns):
    ftype = feature.attrib.get("type", "")
    if ftype in target_types:
        desc = feature.attrib.get("description", "")
        evidence = feature.attrib.get("evidence", "")
        begin = feature.find(".//u:begin", ns)
        end = feature.find(".//u:end", ns)
        begin_pos = begin.attrib.get("position") if begin is not None else None
        end_pos = end.attrib.get("position") if end is not None else None
        
        features.append({
            "type": ftype,
            "description": desc,
            "evidence": evidence,
            "begin": begin_pos,
            "end": end_pos
        })

# 結果を表示
for f in features:
    print(f)
    
output_file =  "domain_and_region_{}.csv".format(uid)
# CSVに書き込み
with open(output_file, mode='w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=["type", "description", "evidence", "begin", "end"])
    writer.writeheader()   # ヘッダー行を書き込む
    writer.writerows(features)

print(f"CSVファイル '{output_file}' に {len(features)} 件のデータを保存しました。")
"""

def main():
    cur = con_db()
    merops_code_mece = create_merops_code_table()
    print("merops_code_mece.")
    #print(merops_code_mece)
    print(len(merops_code_mece))

    #df_domain_region = pd.DataFrame({
    #    'merops_id': [], 'uniprot_id': [], 'p1': [], 'full_aa': [], 'full_aa_len': [], 'type': [], 'description': [], 'evidence': [], 'begin': [], 'end': [], 'region_len': []
    #})
    #domain_number = 0

    df_domain_info_no_exixt = pd.DataFrame({
        'merops_id': [], 'uniprot_id': [], 'p1': [], 'cleave_pattern': [],'full_aa': [], 'full_aa_len': [], 'type': [], 'description': []
    })
    domain_no_exist_number = 0

    num_i = 518#492#0
    for i, merops_id in enumerate(merops_code_mece[num_i:]):
        i = i + num_i 
        df_domain_region = pd.DataFrame({
            'merops_id': [], 'uniprot_id': [], 'p1': [], 'full_aa': [], 'full_aa_len': [], 'type': [], 'description': [], 'evidence': [], 'begin': [], 'end': [], 'region_len': [], 'region_aa': []
        })
        domain_number = 0
    # for i, merops_id in enumerate(['S01.247']):
        print("#"*20 + f"i: {i}_" + f"merops_id: {merops_id}" + "#"*20)
        
        df_cleave_pattern = pd.DataFrame({
            'merops_id': [], 'uniprot_id': [], 'p1': [], 'cleave_pattern': [], 'full_aa': [], 'full_aa_len': []
        })
        cleave_pattern_num = 0
        #display(df_cleave_pattern)
        #print("df_cleave_pattern.")
        #print(df_cleave_pattern)



        print("="*10+"protease: {}".format(merops_id)+"_protease_turn: {}".format(i)+"="*10)

        merops_id = [merops_id]

        cur.execute("SELECT uniprot_acc, p1 FROM cleavage where code=(%s);", merops_id)

        subs = cur.fetchall()
        print("subs: ")
        print(subs)
        print(len(subs))
        len_subs = len(subs)

        num_j = 0 #2060 #2014
        for j, temp_sub in enumerate(subs[num_j:]): # 0, 101, 159. i: 62, j:1468. 
            j = j + num_j
            print("#"*20 + f"i: {i}, " + f"merops_id: {merops_id}" + "#"*20)
            print("#"*30 + f"j: {j}/{len_subs}, " + f"temp_sub: {temp_sub}" + "#"*30)
            uniprot_id = subs[j][0]
            uniprot_id = uniprot_id.strip()
            print(f"uniprot_id: {uniprot_id}")
            content = fr'{uniprot_id}' 

            #pattern = '.*?([A-Z]{1})([A-Z0-9]{2})([A-Z0-9]{3})'
            pattern01 = '.*?([A-Z0-9]{6})$'
            pattern02 = '(^[A-Z0-9]{6}.*?)'
            result01 = re.match(pattern01, content)
            print("result01: ")
            print(result01)
            #result01 = result01[0].strip()
            result02 = re.match(pattern02, content)
            print("result02")
            print(result02)
            #result02 = result02[0].strip()
            if result01 is not None:
                if len(result01[0]) == 6:
                    result = result01
            elif result02 is not None:
                if len(result02[0]) == 6:
                    result = result02
            else:
                ("re error!")
                pass
            #else:
            #    pass
            uniprot_id = result[0]
            print("uniprot_id after re: {}".format(uniprot_id))
            print("result.")
            print(result) 

            if result: #none以外の場合
                #print("result.")
                #print(result) 
                # output:<_sre.SRE_Match object; span=(0, 3), match='hel'>
                #print("result.span().")
                #print(result.span()) 
                # output:(0, 3)
                #print("result.group().")
                #print(result.group()) 
                # output:hel
                pass
            else:   
                #print("result.")
                #print(result) 
                continue
            p1 = subs[j][1]
            print(f"p1: {p1}")

            # UniProt XMLデータのURL
            url = "https://rest.uniprot.org/uniprotkb/{}.xml".format(uniprot_id)
            response = requests.get(url)
            response.raise_for_status()
            # XMLをパース
            root = ET.fromstring(response.content)
            # 名前空間を取得（UniProt XMLは名前空間付き）
            ns = {'u': 'http://uniprot.org/uniprot'}

            # 対象とするtype
            target_types = {"domain", "region of interest", "topological domain", "site", "transmembrane region"}

            # 結果格納用リスト
            features = []

            for feature in root.findall(".//u:feature", ns):
                ftype = feature.attrib.get("type", "")
                if ftype in target_types:
                    desc = feature.attrib.get("description", "")
                    evidence = feature.attrib.get("evidence", "")
                    begin = feature.find(".//u:begin", ns)
                    end = feature.find(".//u:end", ns)
                    begin_pos = begin.attrib.get("position") if begin is not None else None
                    end_pos = end.attrib.get("position") if end is not None else None
                    
                    features.append({
                        "type": ftype,
                        "description": desc,
                        "evidence": evidence,
                        "begin": begin_pos,
                        "end": end_pos
                    })
            print("features: ")
            print(features)
            print(len(features))

            full_aa = aaseq_from_uid(uniprot_id, i, j)
            #full_aa = aaseq_from_uid(uniprot_id, i, j)
            print(f"length of full_aa: "+str(len(full_aa)))
            if len(full_aa) < 3:
                print("full_aa len < 3.")
                continue
            #例外である端っこも取得できるように工夫する
            if p1 - 4 >= 0 and len(full_aa) - p1 >= 4:
                cleave_pattern = full_aa[p1-4:p1+4]
            elif p1 - 4 < 0:
                term = 4 - p1
                cleave_pattern = "-"*term + full_aa[0:8-term]
            elif len(full_aa) - p1 < 4:
                term = 4 - (len(full_aa) - p1)
                cleave_pattern = full_aa[len(full_aa) - (8-term):len(full_aa)] + "-"*term
            else:
                pass            
            # 開裂パターンを表示する
            print("merops_id: {}".format(merops_id))
            print("uniprot_id: {}".format(uniprot_id))
            print("p1: {}".format(p1))
            print("cleave_pattern: {}".format(cleave_pattern))
            print("full_aa: {}".format(full_aa))
            
            df_cleave_pattern.loc[f'{cleave_pattern_num}'] = [";".join(merops_id), uniprot_id, p1, cleave_pattern, full_aa, len(full_aa)]
            filename_cp = f'./proteases/cleave_pattern_one_letter_aa_{merops_id[0]}_withDomain.csv'
            df_cleave_pattern.to_csv(filename_cp)
            #print("df_cleave_pattern.")
            #print(df_cleave_pattern)
            cleave_pattern_num += 1


            print("features is None?")    
            #if features is None:
            if len(features) == 0:
                print("features is None!")
                df_domain_info_no_exixt.loc[f'{domain_no_exist_number}'] = [";".join(merops_id), uniprot_id, p1, cleave_pattern, full_aa, len(full_aa), None, None]
                filename_domain_no_exist = f'./domains/domain_no_exist.csv'
                df_domain_info_no_exixt.to_csv(filename_domain_no_exist)
                domain_no_exist_number += 1
                #'merops_id': [], 'uniprot_id': [], 'p1': [], 'full_aa': []
                continue 

            print("features is Exist.")    
            for f in features:
                print("p1: {}".format(p1))
                print(" f[begin]: {}".format(f["begin"]))
                print(" f[end]: {}".format(f["end"]))
                print("f: ")
                print(f)
                if f["begin"] is None or f["end"] is None:
                    print("f[begin] or f[end] is None!")
                    df_domain_info_no_exixt.loc[f'{domain_no_exist_number}'] = [";".join(merops_id), uniprot_id, p1, cleave_pattern, full_aa, len(full_aa), f["type"], f["description"]]
                    filename_domain_no_exist = f'./domains/domain_no_exist.csv'
                    df_domain_info_no_exixt.to_csv(filename_domain_no_exist)
                    domain_no_exist_number += 1
                    #'merops_id': [], 'uniprot_id': [], 'p1': [], 'full_aa': []
                    break
                if p1 >= int(f["begin"]) or p1 <= int(f["end"]):
                    print(" f[begin]: {}".format(f["begin"]))
                    print("p1: {}".format(p1))
                    print("f: ")
                    print(f)
                    region_len = int(f["end"]) - int(f["begin"]) + 1
                    df_domain_region.loc[f'{domain_number}'] = [";".join(merops_id), uniprot_id, p1, full_aa, len(full_aa), f["type"], f["description"], f["evidence"], f["begin"], f["end"], region_len, full_aa[int(f["begin"])-1: int(f["end"])-1]]
                    filename_dr = f'./domains/domain_and_region_{merops_id[0]}.csv'
                    df_domain_region.to_csv(filename_dr)
                    # 'type': [], 'description': [], 'evidence': [], 'begin': [], 'end': []
                    domain_number += 1


            
            #output_file =  "./domains/domain_and_region_{}_{}.csv".format(merops_id[0], uniprot_id)
            #with open(output_file, mode='w', newline='', encoding='utf-8') as f:
            #    writer = csv.DictWriter(f, fieldnames=["type", "description", "evidence", "begin", "end"])
            #    writer.writeheader()   # ヘッダー行を書き込む
            #    writer.writerows(features)
            

            filename_domain_no_exist = './domains/domain_no_exist.csv'
            df_domain_info_no_exixt.to_csv(filename_domain_no_exist)
            

        filename_cp = f'./proteases/cleave_pattern_one_letter_aa_{merops_id[0]}_withDomain.csv'
        df_cleave_pattern.to_csv(filename_cp)
        print("df_cleave_pattern.")
        print(df_cleave_pattern)   

        filename_dr = f'./domains/domain_and_region_{merops_id[0]}.csv'
        df_domain_region.to_csv(filename_dr)
        print("df_domain_region.")
        print(df_domain_region)
        

    print("END.")
    sys.exit()



    print("main END.")
    
def con_db():
    # コネクションの作成
    conn = mydb.connect(
        host='localhost',
        port=3306,
        user='foo',
        password='bar',
        database='meropsweb121'
    )
    # DB操作用にカーソルを作成
    cur = conn.cursor()
    return cur

def safe_read_response(response, retries=3, delay=2):
    """IncompleteReadが出てもリトライして読み取る"""
    for attempt in range(retries):
        try:
            return response.read()
        except http.client.IncompleteRead as e:
            print(f"IncompleteRead error: {e}. Retrying... ({attempt+1}/{retries})")
            if attempt < retries - 1:
                time.sleep(delay)
            else:
                raise  # リトライ失敗時は再度例外送出


def aaseq_from_uid(uid, protease_turn, substrate_turn):
    df = pd.DataFrame(np.arange(3).reshape(1, 3), columns=['uniprotKB_accession', 'function', 'sequence'], index=['protease'+str(protease_turn)+'_substrate'+str(substrate_turn)])

    for column_name in df:
        df[column_name] = df[column_name].astype(str)

    df['uniprotKB_accession'][0] = uid

    #display(df)    

    url = "https://www.uniprot.org/uniprot/" + uid + ".xml"
    f = urlopen(url)
    # xml = f.read()
    xml = safe_read_response(f)
    root = etree.fromstring(xml)
    
    #以下のコードは下の説明を参照
    function = root.find('./entry/comment[@type="function"]', root.nsmap)
    if function==None:
        print("function was not detected.")
        pass
    else:
        df["function"][0] = function[0].text
        #print(function[0].text+"¥n")

    sequence = root.find('./entry/sequence', root.nsmap) 
    if sequence==None: 
        print("sequence was not detected.")
        pass 
    else: 
        df["sequence"][0] = sequence.text 
        #print(sequence.text+"¥n")

    #display(df) 
    #df0=pd.concat([df0, df], axis=0)
    #display(df0)
    print(df)
    print(df["sequence"][0]) 

    return df["sequence"][0]



def create_merops_code_table():
    csv_file = open("./testdata/merops_code_mece.csv", "r", encoding="ms932", errors="", newline="" )
    #リスト形式
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

    #print(f)
    #print(len(f))

    i = 0
    for row in f:
        print(row)
        merops_code_mece = row
        i = i + 1
        if i == 1:
            break

    #print("merops_code_mece")
    #print(merops_code_mece)
    #len(merops_code_mece)
    return merops_code_mece

if __name__ == '__main__':
    start = time.time()
    
    main()

    end = time.time()
    time_diff = end - start
    print(time_diff)
    print(time_diff/60)
    print(time_diff/60/60)
    print(time_diff/60/60/24)
    
    dt_now = datetime.datetime.now()
    print(dt_now)
    print("END.")

