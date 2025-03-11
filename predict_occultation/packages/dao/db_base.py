import collections
import configparser
import os
import warnings

from sqlalchemy import MetaData, Table, create_engine
from sqlalchemy import exc as sa_exc
from sqlalchemy.pool import NullPool
from sqlalchemy.sql import and_, select


class MissingDBURIException(Exception):
    pass


class DBBase:

    con = None

    def get_db_uri(self):

        try:
            db_uri = os.environ["DB_ADMIN_URI"]
            return db_uri
        except:
            raise MissingDBURIException(
                "Required environment variable with URI to access the database."
                "example DB_ADMIN_URI=postgresql+psycopg2://USER:PASS@HOST:PORT/DB_NAME"
            )

    def get_db_engine(self):
        engine = create_engine(self.get_db_uri(), poolclass=NullPool)

        return engine

    def get_con(self):
        if self.con is None:
            engine = self.get_db_engine()
            self.con = engine.connect()

        return self.con

    def get_table(self, tablename, schema=None):

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=sa_exc.SAWarning)

            engine = self.get_db_engine()
            tbl = Table(tablename, MetaData(engine), autoload=True, schema=schema)
            return tbl

    def fetch_all_dict(self, stm):

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=sa_exc.SAWarning)

            engine = self.get_db_engine()
            with engine.connect() as con:

                queryset = con.execute(stm)

                rows = []
                for row in queryset:
                    d = dict(collections.OrderedDict(row))
                    rows.append(d)

                return rows

    def fetch_one_dict(self, stm):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=sa_exc.SAWarning)

            engine = self.get_db_engine()
            with engine.connect() as con:

                queryset = con.execute(stm).fetchone()

                if queryset is not None:
                    d = dict(collections.OrderedDict(queryset))
                    return d
                else:
                    return None

    def fetch_scalar(self, stm):
        engine = self.get_db_engine()
        with engine.connect() as con:
            return con.execute(stm).scalar()

    # def get_job_by_id(self, id):

    #     tbl = self.get_table(tablename="des_astrometryjob")
    #     stm = select(tbl.c).where(and_(tbl.c.id == int(id)))

    #     return self.fetch_one_dict(stm)

    def import_with_copy_expert(self, sql, data):
        """
            This method is recommended for importing large volumes of data. using the postgresql COPY method.

            The method is useful to handle all the parameters that PostgreSQL makes available
            in COPY statement: https://www.postgresql.org/docs/current/sql-copy.html

            it is necessary that the from clause is reading from STDIN.

            example:
            sql = COPY <table> (<columns) FROM STDIN with (FORMAT CSV, DELIMITER '|', HEADER);

            Parameters:
                sql (str): The sql statement should be in the form COPY table '.
                data (file-like ): a file-like object to read or write
            Returns:
                rowcount (int):  the number of rows that the last execute*() produced (for DQL statements like SELECT) or affected (for DML statements like UPDATE or INSERT)

        References:
            https://www.psycopg.org/docs/cursor.html#cursor.copy_from
            https://stackoverflow.com/questions/30050097/copy-data-from-csv-to-postgresql-using-python
            https://stackoverflow.com/questions/13125236/sqlalchemy-psycopg2-and-postgresql-copy
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=sa_exc.SAWarning)

            connection = self.get_db_engine().raw_connection()
            try:
                cursor = connection.cursor()
                cursor.copy_expert(sql, data)
                connection.commit()

                cursor.close()
                return cursor.rowcount
            except Exception as e:
                connection.rollback()
                raise (e)
            finally:
                connection.close()
