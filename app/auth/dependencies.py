from fastapi import Depends, HTTPException, status
from fastapi.security import OAuth2PasswordBearer
from jose import JWTError, jwt
from datetime import datetime, timedelta
from typing import Optional

from sqlalchemy.orm import Session

from ..db.session import get_db
from ..models.Users import Users, UserDB
from .jwt_handler import SECRET_KEY, ALGORITHM, TokenData

oauth2_scheme = OAuth2PasswordBearer(tokenUrl="auth/token")

def get_current_user(token: str = Depends(oauth2_scheme), db = Depends(get_db)):
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        email: str = payload.get("sub")
        if email is None:
            raise credentials_exception
        token_data = TokenData(email=email)
    except JWTError:
        raise credentials_exception
    
    user = db.query(Users).filter(Users.email == email).first()
    if user is None:
        raise credentials_exception
    return user

def get_current_active_user(current_user: Users = Depends(get_current_user)):
    if not current_user.is_active:
        raise HTTPException(status_code=400, detail="Inactive user")
    return current_user

def admin_required(current_user: Users = Depends(get_current_active_user)):
    if current_user.role != "admin":
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Admin privileges required"
        )
    return current_user

def get_password_hash(password: str) -> str:
    from app.auth.jwt_handler import password_hash
    return password_hash.hash(password)

def verify_password(plain_password: str, hashed_password: str) -> bool:
    from app.auth.jwt_handler import password_hash
    return password_hash.verify(plain_password, hashed_password)

def authenticate_user(db_session: Session, email: str, password: str):
    user = db_session.query(Users).filter(Users.email == email).first()
    if not user:
        return False
    
    user_db = db_session.query(UserDB).filter(UserDB.user_id == user.id).first()
    if not user_db:
        return False
        
    if not verify_password(password, user_db.hashed_password):
        return False
        
    return user_db