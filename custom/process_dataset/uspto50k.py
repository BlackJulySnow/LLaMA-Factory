import lmdb

env = lmdb.open('custom\dataset\USPTO50K\train.lmdb')
with env.begin() as txn:
    cursor = txn.cursor()
    for key, value in cursor:
        print(key, value)