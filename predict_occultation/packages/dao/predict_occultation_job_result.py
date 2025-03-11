from datetime import timedelta

from dao.db_base import DBBase
from sqlalchemy.sql import and_, delete, select
from sqlalchemy import func
from typing import List

class PredictOccultationJobResultDao(DBBase):
    def __init__(self):
        super(PredictOccultationJobResultDao, self).__init__()

        self.tbl = self.get_table("tno_predictionjobresult")

    def delete_by_job_id(self, job_id):

        stm = delete(self.tbl).where(and_(self.tbl.c.job_id == job_id))

        engine = self.get_db_engine()
        with engine.connect() as con:
            rows = con.execute(stm)

            return rows

    def insert(self, data):
        engine = self.get_db_engine()
        with engine.connect() as con:
            result = con.execute(self.tbl.insert(), data)
            return result.inserted_primary_key[0]

    def update(self, id, data):

        if "exec_time" in data:
            data["exec_time"] = timedelta(seconds=data["exec_time"])
        if "pre_occ_exec_time" in data:
            data["pre_occ_exec_time"] = timedelta(seconds=data["pre_occ_exec_time"])
        if "calc_path_coeff_exec_time" in data:
            data["calc_path_coeff_exec_time"] = timedelta(
                seconds=data["calc_path_coeff_exec_time"]
            )
        if "ing_occ_exec_time" in data:
            data["ing_occ_exec_time"] = timedelta(seconds=data["ing_occ_exec_time"])

        engine = self.get_db_engine()
        with engine.connect() as con:
            result = con.execute(self.tbl.update().where(self.tbl.c.id == id), data)
            return result.rowcount

    def by_job_id(self, job_id):
        stm = select(self.tbl.c).where(and_(self.tbl.c.job_id == job_id))

        return self.fetch_all_dict(stm)

    def count_by_job_id(self, job_id, status: List = None):
        stm = select([func.count()]).select_from(self.tbl).where(and_(self.tbl.c.job_id == job_id))
        if status:
            stm = stm.where(self.tbl.c.status.in_(status))

        return self.fetch_scalar(stm)