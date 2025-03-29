use solana_program::{
    account_info::AccountInfo,
    entrypoint::ProgramResult,
    pubkey::Pubkey,
};
use super::*;

#[test]
fn test_create_data_ownership() {
    let program_id = Pubkey::new_unique();
    let mut data = vec![0; 64];
    let mut lamports = 0;
    let owner = Pubkey::new_unique();
    let data_hash = [1; 32];
    
    let data_account = AccountInfo::new(
        &Pubkey::new_unique(),
        &mut data,
        &mut lamports,
        &mut vec![],
        &program_id,
        false,
        false,
        false,
    );

    let accounts = vec![data_account];
    
    let mut instruction_data = vec![0];
    instruction_data.extend_from_slice(&owner.to_bytes());
    instruction_data.extend_from_slice(&data_hash);
    
    let result = process_instruction(&program_id, &accounts, &instruction_data);
    assert!(result.is_ok());
    
    assert_eq!(&data[0..32], &owner.to_bytes());
    assert_eq!(&data[32..64], &data_hash);
}

#[test]
fn test_transfer_ownership() {
    let program_id = Pubkey::new_unique();
    let mut data = vec![0; 64];
    let mut lamports = 0;
    let old_owner = Pubkey::new_unique();
    let new_owner = Pubkey::new_unique();
    
    // Initialize data with old owner
    data[0..32].copy_from_slice(&old_owner.to_bytes());
    
    let data_account = AccountInfo::new(
        &Pubkey::new_unique(),
        &mut data,
        &mut lamports,
        &mut vec![],
        &program_id,
        false,
        false,
        false,
    );
    
    let owner_account = AccountInfo::new(
        &old_owner,
        &mut vec![],
        &mut lamports,
        &mut vec![],
        &program_id,
        true,
        true,
        false,
    );
    
    let accounts = vec![data_account, owner_account];
    
    let mut instruction_data = vec![1];
    instruction_data.extend_from_slice(&new_owner.to_bytes());
    
    let result = process_instruction(&program_id, &accounts, &instruction_data);
    assert!(result.is_ok());
    
    assert_eq!(&data[0..32], &new_owner.to_bytes());
} 