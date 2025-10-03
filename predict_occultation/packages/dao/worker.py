import datetime
import traceback

from dao.db_base import DBBase
from sqlalchemy import update
import datetime

class WorkerDao(DBBase):
    def __init__(self):
        super(WorkerDao, self).__init__()
        self.tbl = self.get_table("predict_occultation_workersheartbeat")

    def initialize_heartbeat(self, worker_name):
        with self.get_db_engine().connect() as db:
            try:
                heartbeat = db.execute(self.tbl.select().where(self.tbl.c.worker == worker_name)).first()
                if not heartbeat:
                    stmt = self.tbl.insert().values(
                        worker=worker_name,
                        started_at=datetime.datetime.now(tz=datetime.timezone.utc),
                        updated_at=datetime.datetime.now(tz=datetime.timezone.utc),
                        uptime=0,
                    )
                    db.execute(stmt)
                    db.commit()
                    # print(f"Heartbeat initialized for worker {worker_name}.")
                else:
                    stmt = (
                        update(self.tbl)
                        .where(self.tbl.c.worker == heartbeat.worker)
                        .values(
                            started_at = datetime.datetime.now(tz=datetime.timezone.utc),
                            updated_at=datetime.datetime.now(tz=datetime.timezone.utc),
                            uptime=0
                        )
                    )
                    db.execute(stmt)
                    db.commit()
                    # self.log.debug(f"Heartbeat updated for worker {worker_name}.")
            except Exception as e:
                print(f"Error initializing/updating heartbeat: {e}")
                # add traceback for debugging
                print(traceback.format_exc())
                db.rollback()
            finally:
                db.close()



    def send_heartbeat(self, worker_name):
        with self.get_db_engine().connect() as db:
            try:
                heartbeat = db.execute(self.tbl.select().where(self.tbl.c.worker == worker_name)).first()
                if not heartbeat:
                    self.initialize_heartbeat(worker_name)
                    return
                stmt = (
                    update(self.tbl)
                    .where(self.tbl.c.worker == worker_name)
                    .values(
                        uptime=int((
                            datetime.datetime.now(tz=datetime.timezone.utc)
                            - heartbeat.started_at
                        ).total_seconds()),
                        updated_at=datetime.datetime.now(tz=datetime.timezone.utc)
                    )
                )

                db.execute(stmt)
                db.commit()
                # print(f"Heartbeat sent for worker {worker_name}.")
                return True
            except Exception as e:
                msg = f"Error sending heartbeat for worker {worker_name}: {e}"
                db.rollback()
                raise Exception(msg)

