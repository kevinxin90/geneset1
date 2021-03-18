import requests
import csv
import os

TYPE_TO_MYGENE_FIELD_MAPPING = {
    "pathway": "pathway",
    "bp": "go.BP",
    "mf": "go.MF",
    "cc": "go.CC"
}

pathway_ids = set()

def count_pathway_participants(rec_id):
    [db, record_id] = rec_id.split(':')
    url = 'https://mygene.info/v3/query?q=pathway.' + db + '.id:"' + record_id + '"&size=0'
    doc = requests.get(url).json()
    return doc.get("total", 0)

def count_go_participants(_type, record_id):
    url = 'https://mygene.info/v3/query?q=go.' + _type.upper() + '.id:GO\:' + record_id.split(':')[-1]  + '&size=0'
    print(url)
    doc = requests.get(url).json()
    return doc.get("total", 0)

def parse_pathway_response(doc):
    res = []
    for _doc in doc['hits']:
        for db, info in _doc['pathway'].items():
            if isinstance(info, dict):
                info = [info]
            for record in info:
                _id = db + ':' + str(record['id'])
                if _id not in pathway_ids:
                    pathway_ids.add(_id)
                    res.append(
                        {
                            '_id': _id,
                            db: record['id'],
                            'name': record['name'],
                            'type': 'pathway'
                        }
                    )
    return res

def parse_go_respnse(doc, _type):
    res = []
    for _doc in doc['hits']:
        info = _doc['go'][_type.upper()]
        if isinstance(info, dict):
            info = [info]
        for record in info:
            _id = record['id']
            if _id not in pathway_ids:
                pathway_ids.add(_id)
                yield {'_id': _id,
                        'go': _id,
                        'name': record['term'],
                        'type': _type}

def load_from_mygene(semanticType):
    mygene_field = TYPE_TO_MYGENE_FIELD_MAPPING[semanticType]
    url = 'http://mygene.info/v3/query?q=_exists_:' + mygene_field + '&fields=' + mygene_field + '&fetch_all=TRUE'
    cnt = 0
    total = 1
    results = []
    while cnt < total:
        doc = requests.get(url).json()
        if total == 1:
            total = doc['total']
        cnt += len(doc['hits'])
        url = 'http://mygene.info/v3/query?scroll_id=' + doc['_scroll_id']
        if semanticType == "pathway":
            results += parse_pathway_response(doc)
        else:
            results += parse_go_respnse(doc, semanticType)
    if semanticType == "pathway":
        for rec in results:
            rec['num_of_participants'] = count_pathway_participants(rec['_id'])
            yield rec
    else:
        for rec in results:
            rec['num_of_participants'] = count_go_participants(semanticType, rec['_id'])
            yield rec

def load_data(data_folder):
    # load cellular component and biological proces info from semmed neo4j
    nodes_path = os.path.join(data_folder, "nodes_neo4j.csv")
    with open(nodes_path) as f:
        csv_reader = csv.reader(f, delimiter=',')
        next(csv_reader)
        for _item in csv_reader:
            semantic = _item[-2]
            if semantic == 'biological_process_or_activity':
                yield {'_id': _item[-1],
                       'umls': _item[-1].split(':')[-1],
                       'type': 'bp',
                       'name': _item[1]}
            elif semantic == 'cell_component':
                yield {'_id': _item[-1],
                       'umls': _item[-1].split(':')[-1],
                       'type': 'cc',
                       'name': _item[1]}
    # load pathways
    for rec in load_from_mygene("pathway"):
        yield rec
    print("pathway done")
    # load biological process
    for rec in load_from_mygene("bp"):
        yield rec
    print("bp done")
    # load molecular function
    for rec in load_from_mygene("mf"):
        yield rec
    print("mf done")
    # load cellular component
    for rec in load_from_mygene("cc"):
        yield rec
    print("cc done")