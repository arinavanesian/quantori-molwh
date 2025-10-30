from jose import JWTError, jwt
from datetime import datetime, timezone, timedelta
from fastapi import HTTPException, status, Depends
from pwdlib import PasswordHash
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm

from pydantic import BaseModel
from models import Users, UserDB
from db import engine, SessionLocal

# openssl rand -hex 32
SECRET_KEY = "bc0db4a2d1062ecaff76294bf8609fa33c4c401bac322e5bcbd07c5c5e048885"
ALGORITHM = "HS256"

password_hash = PasswordHash.recommended()
oauth_scheme = OAuth2PasswordBearer(tokenUrl="token")
class TokenData(BaseModel):
    email:str | None = None
    role:str | None = None




def create_token(data:dict, SECRET_KEY:str, ALGORITHM:str, expires_delta: timedelta | None = None,):
    to_encode  = data.copy()
    expire = datetime.now(timezone.utc) + (expires_delta or timedelta(minutes=30))
    to_encode.update({"exp":expire})
    encode_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encode_jwt


def verify_token(token:str, SECRET_KEY:str, ALGORITHM:str):
    credential_exception = HTTPException(status_code =status.HTTP_401_UNAUTHORIZED,
                                         detail = "Token invalid! Authorization not passed!",
                                         headers = {"WWW-Authenticate":"Bearer"},)
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        email:str = payload.get("email")
        role:str = payload.get("role")
        if not email or not role:
            raise credential_exception
        token_data = TokenData(email=email, role=role)
        return token_data
    except jwt.ExpiredSignatureError:
        raise JWTError("Token has expired")
    except jwt.InvalidTokenError:
        raise JWTError("Invalid token")
    

def is_admin(token_data:TokenData = Depends(verify_token)):
    try:
        if token_data.role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN,
                                detail="Admin privileges required!")
    except JWTError:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN,
                            detail="Could not validate credentials")
    return token_data
