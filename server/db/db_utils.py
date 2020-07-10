import typing

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from server.db.cellxgene_orm import Base


class DbUtils:
    def __init__(self, database_uri: str = "postgresql://postgres:test_pw@localhost:5432"):
        self.session = DBSessionMaker(database_uri).session()
        self.engine = self.session.get_bind()

    def get(self, table: Base, entity_id: typing.Union[str, typing.Tuple[str]]) -> typing.Union[Base, None]:
        """
        Query a table row by its primary key
        :param table: SQLAlchemy Table to query
        :param entity_id: Primary key of desired row
        :return: SQLAlchemy Table object, None if not found
        """
        return self.session.query(table).get(entity_id)

    def query(self, table_args: typing.List[Base], filter_args: typing.List[bool] = None) -> typing.List[Base]:
        """
        Query the database using the current DB session
        :param table_args: List of SQLAlchemy Tables to query/join
        :param filter_args: List of SQLAlchemy filter conditions
        :return: List of SQLAlchemy query response objects
        """
        return (
            self.session.query(*table_args).filter(*filter_args).all()
            if filter_args
            else self.session.query(*table_args).all()
        )


class DBSessionMaker:
    def __init__(self, database_uri):
        self.engine = create_engine(database_uri, connect_args={"connect_timeout": 5})
        self.session_maker = sessionmaker(bind=self.engine)

    def session(self, **kwargs):
        return self.session_maker(**kwargs)
