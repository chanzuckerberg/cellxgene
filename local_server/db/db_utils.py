import typing
import uuid

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from local_server.db.cellxgene_orm import Base, CellxGeneDataset, CellxGeneUser


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

    def query_for_most_recent(self, table: Base, filter_args: typing.List[bool] = None) -> Base:
        try:
            return self.session.query(table).filter(*filter_args).order_by(table.created_at.desc()).limit(1).all()[0]
        except IndexError:
            return None

    def get_or_create_dataset(self, dataset_name):
        try:
            dataset_id = self.query(table_args=[CellxGeneDataset], filter_args=[CellxGeneDataset.name == dataset_name])[
                0
            ].id
        except IndexError:
            dataset_id = uuid.uuid4()
            dataset = CellxGeneDataset(id=dataset_id, name=dataset_name)
            self.session.add(dataset)
            self.session.commit()
        return str(dataset_id)

    def get_or_create_user(self, user_id):
        try:
            user_id = self.query(table_args=[CellxGeneUser], filter_args=[CellxGeneUser.id == user_id])[0].id
        except IndexError:
            user = CellxGeneUser(id=user_id)
            self.session.add(user)
            self.session.commit()
        return str(user_id)


class DBSessionMaker:
    def __init__(self, database_uri):
        self.engine = create_engine(database_uri, connect_args={"connect_timeout": 5})
        self.session_maker = sessionmaker(bind=self.engine)

    def session(self, **kwargs):
        return self.session_maker(**kwargs)
