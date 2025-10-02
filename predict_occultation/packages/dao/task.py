import datetime
from dao.db_base import DBBase
from sqlalchemy import update
from sqlalchemy.sql import and_,  select
from sqlalchemy.exc import OperationalError
from sqlalchemy.dialects import postgresql
from enum import StrEnum

class PredictionState(StrEnum):
    PENDING = "PENDING"
    PREPARING = "PREPARING"
    READY_FOR_RUN = "READY_FOR_RUN"
    SUBMITTING = "SUBMITTING"
    QUEUED = "QUEUED"
    RUNNING = "RUNNING"
    WAITING_RESULTS = "WAITING_RESULTS"
    INGESTING = "INGESTING"
    DONE = "DONE"
    FAILED = "FAILED"
    STALLED = "STALLED"
    ABORTED = "ABORTED"

class TaskDao(DBBase):
    def __init__(self):
        super(TaskDao, self).__init__()
        self.base_delay = 60  # Base delay in seconds for retries
        self.tbl = self.get_table("predict_occultation_predictiontask")

    def debug_query(self, stm, with_parameters=False):
        sql = str(
            stm.compile(
                dialect=postgresql.dialect(),
                compile_kwargs={"literal_binds": with_parameters},
            )
        )

        # Remove new lines
        sql = sql.replace("\n", " ").replace("\r", "")
        return sql


    def get_next_task(self, db_session, state_to_process):
        try:
            stm = (
                select(self.tbl.c)
                .where(
                    self.tbl.c.state == state_to_process,
                    self.tbl.c.aborted == False,
                    (self.tbl.c.next_retry_at == None) | (self.tbl.c.next_retry_at <= datetime.datetime.now(tz=datetime.timezone.utc))
                )
                .order_by(self.tbl.c.priority.desc(), self.tbl.c.created_at.asc())
                .with_for_update(skip_locked=True)
                .limit(1)
            )

            # print(f"SQL: {self.debug_query(stm, with_parameters=True)}")

            result = db_session.execute(stm).first()
            return result
    
        except OperationalError as e:
            print(f"OperationalError: {e}")
            # self.log.warning(f"Erro operacional ao buscar task (provavelmente bloqueio): {e}")
            db_session.rollback()
            return None
        except Exception as e:
            print(f"Error fetching next task: {e}")
            # self.log.error(f"Erro ao buscar a próxima task: {e}")
            db_session.rollback()
            return None
        
    def update_task_status(self, db_session, task_id:int, new_state:str):
        try:
            stmt = (
                update(self.tbl)
                .where(self.tbl.c.id == task_id)
                .values(state=new_state, updated_at=datetime.datetime.now(tz=datetime.timezone.utc))
            )

            db_session.execute(stmt)
            db_session.commit()
            print(f"Changed task state to: {new_state}")
            return True
        except Exception as e:
            msg = f"Erro ao atualizar o status da task {task_id}: {e}"
            db_session.rollback()
            raise Exception(msg)
        
    def get_task_by_id(self, db_session, task_id:int):
        try:
            stm = select(self.tbl.c).where(self.tbl.c.id == task_id)
            result = db_session.execute(stm).first()
            return result
        except Exception as e:
            msg = f"Erro ao buscar a task {task_id}: {e}"
            db_session.rollback()
            raise Exception(msg)
        
    def mark_task_failed(self, db_session, task_id, error_message):
        try:
            task = self.get_task_by_id(db_session, task_id)

            attempt_count = task.attempt_count + 1
            if attempt_count >= task.max_retries:
                state = PredictionState.STALLED
                next_retry_at = None
                print(f"Task {task.id} marcada como STALLED após atingir o número máximo de retries.")
            else:
                # Exponential backoff
                delay_seconds = self.base_delay * (2 ** (task.attempt_count - 1)) 
                state = PredictionState.PENDING
                next_retry_at = datetime.datetime.now(tz=datetime.timezone.utc) + datetime.timedelta(seconds=delay_seconds)
                print(f"Task {task.id} marcada como FAILED. Próximo retry em {delay_seconds} segundos.")

            updated_at = datetime.datetime.now(tz=datetime.timezone.utc)
            stmt = (
                update(self.tbl)
                .where(self.tbl.c.id == task_id)
                .values(
                    state=state,
                    last_error=error_message,
                    attempt_count=attempt_count,
                    next_retry_at=next_retry_at,
                    updated_at=updated_at
                )
            )
            db_session.execute(stmt)
            db_session.commit()
        except Exception as e:
            db_session.rollback()
            msg = f"Erro ao marcar a task {task_id} como FAILED: {e}"
            raise Exception(msg)

