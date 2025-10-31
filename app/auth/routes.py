from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from datetime import timedelta
from sqlalchemy.orm import Session

from ..db.session import get_db
from .dependencies import (
    authenticate_user,
    get_password_hash, get_current_active_user)
from .jwt_handler import create_access_token, ACCESS_TOKEN_EXPIRE_MINUTES
from ..models.Users import Users, UserDB

router = APIRouter(prefix="/auth", tags=["authentication"])

@router.post("/token")
async def login_for_access_token(
    form_data: OAuth2PasswordRequestForm = Depends(),
    db: Session = Depends(get_db)
):
    user = authenticate_user(db, form_data.username, form_data.password)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect email or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
    access_token_expires = timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    access_token = create_access_token(
        data={"sub": user.user.email}, expires_delta=access_token_expires
    )
    
    return {
        "access_token": access_token,
        "token_type": "bearer",
        "user_email": user.user.email,
        "user_role": user.user.role
    }

@router.post("/register")
async def register_user(
    email: str,
    password: str,
    role: str = "user",
    db: Session = Depends(get_db)
):
    existing_user = db.query(Users).filter(Users.email == email).first()
    if existing_user:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Email already registered"
        )
    
    new_user = Users(email=email, role=role)
    db.add(new_user)
    db.flush()
    
    hashed_password = get_password_hash(password)
    new_user_db = UserDB(user_id=new_user.id, hashed_password=hashed_password)
    db.add(new_user_db)
    
    db.commit()
    
    return {"message": "User created successfully", "email": email}

@router.get("/me")
async def read_users_me(current_user: Users = Depends(get_current_active_user)):
    return {
        "email": current_user.email,
        "role": current_user.role,
        "is_active": current_user.is_active
    }