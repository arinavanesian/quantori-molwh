from sqlalchemy.orm import sessionmaker, Session
from contextlib import contextmanager


from .engine import engine
from typing import Generator
SessionLocal = sessionmaker(bind = engine, autocommit = False, autoflush = False, expire_on_commit=False)

# @contextmanager
def get_db()->Generator[Session, None, None]:
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
