from sqlalchemy.orm import sessionmaker, Session
from contextlib import contextmanager

# from .. import db
# from app.db import engine
from .engine import engine
# from app.db.engine import engine

SessionLocal = sessionmaker(bind = engine, autocommit = False, autoflush = False, expire_on_commit=False)

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

@contextmanager
def session_scope():
    session = SessionLocal()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()
