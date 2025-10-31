from app.db.Base import Base
from app.db.engine import engine
from app.db.session import session_scope
from sqlalchemy.orm import Mapped, mapped_column, relationship
from sqlalchemy import (
    String, Integer, Column, 
    DateTime, ForeignKey)
from sqlalchemy import Enum as SQLEnum

from enum import Enum
from datetime import datetime

from typing import List

import uuid as uuid

class UserRole(str, Enum):
    ADMIN = "admin"
    USER = "user"
    
class Users(Base):
    
    """one-to-one relationship with the UserDB
    one-to-many relationship with the Query"""
    __tablename__ = "users"
    id : Mapped[int] = mapped_column(Integer, primary_key=True, index=True)
    email: Mapped[str] = mapped_column(String, unique=True, index=True, nullable=False)
    role: Mapped[UserRole] = mapped_column(SQLEnum(UserRole), default=UserRole.USER, nullable=False )
    is_active: Mapped[bool] = mapped_column(default=True, nullable=False)

    user_db: Mapped["UserDB"] = relationship("UserDB", back_populates="user", uselist=False)
    queries: Mapped[List["Query"]] = relationship("Query", back_populates="user")


    def __repr__(self):
        return f"<User(email={self.email}, role={self.role})>"

class UserDB(Base):
    """Stores the users hashed password with a
    one-to-one relationship with the Users"""
    __tablename__ = "user_db"
    user_id: Mapped[int] = mapped_column(ForeignKey('users.id'), primary_key=True, index = True, nullable=False)
    hashed_password: Mapped[str] = mapped_column(String, nullable=False)

    user: Mapped["Users"] = relationship("Users", back_populates="user_db")

class MoleculeWarehouse(Base):
    __tablename__ = "molecule_store"
    
    uuid: Mapped[str] = mapped_column(String, primary_key=True, 
                                      index=True, nullable=False,
                                      default = lambda: str(uuid.uuid4()))
    smiles: Mapped[str] = mapped_column(String, unique=True, index=True, nullable=False)
    iupac_name: Mapped[str] = mapped_column(String, nullable=True)

    queries: Mapped[List["Query"]] = relationship(
        "Query", 
        back_populates="molecule",
        cascade="all, delete-orphan"
    )
    
class Query(Base):
    """
    A table tracking the User querries, later to be
    used for analyzing the scientific interests of users
     many-to-one relationship with the Users"""
    __tablename__ = "query"
    query_id:Mapped[int] = mapped_column(Integer, primary_key=True, index=True)
    user_id:Mapped[int]= mapped_column(ForeignKey("users.id"), nullable=False)
    molecule_smiles:Mapped[str]= mapped_column(ForeignKey("molecule_store.uuid"), nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime, default= datetime.now, nullable=False)

    user: Mapped["Users"] = relationship("Users", back_populates="queries")
    molecule: Mapped["MoleculeWarehouse"] = relationship("MoleculeWarehouse", back_populates="queries")

