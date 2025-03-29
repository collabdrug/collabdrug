use solana_program::{
    account_info::{next_account_info, AccountInfo},
    entrypoint,
    entrypoint::ProgramResult,
    pubkey::Pubkey,
};

// Declare the program's entrypoint
entrypoint!(process_instruction);

// Program entrypoint implementation
pub fn process_instruction(
    program_id: &Pubkey,
    accounts: &[AccountInfo],
    instruction_data: &[u8],
) -> ProgramResult {
    // Parse the instruction data
    let instruction = DataOwnershipInstruction::unpack(instruction_data)?;

    match instruction {
        DataOwnershipInstruction::CreateDataOwnership { owner, data_hash } => {
            process_create_data_ownership(program_id, accounts, owner, data_hash)
        }
        DataOwnershipInstruction::TransferOwnership { new_owner } => {
            process_transfer_ownership(program_id, accounts, new_owner)
        }
    }
}

// Define the program's instructions
#[derive(Debug)]
enum DataOwnershipInstruction {
    CreateDataOwnership {
        owner: Pubkey,
        data_hash: [u8; 32],
    },
    TransferOwnership {
        new_owner: Pubkey,
    },
}

impl DataOwnershipInstruction {
    fn unpack(input: &[u8]) -> Result<Self, ProgramError> {
        let (&tag, rest) = input.split_first().ok_or(ProgramError::InvalidInstructionData)?;
        Ok(match tag {
            0 => {
                let (owner, rest) = rest.split_at(32);
                let owner = owner.try_into().map_err(|_| ProgramError::InvalidInstructionData)?;
                let data_hash = rest.try_into().map_err(|_| ProgramError::InvalidInstructionData)?;
                Self::CreateDataOwnership {
                    owner: Pubkey::new_from_array(owner),
                    data_hash,
                }
            }
            1 => {
                let (new_owner, _) = rest.split_at(32);
                let new_owner = new_owner.try_into().map_err(|_| ProgramError::InvalidInstructionData)?;
                Self::TransferOwnership {
                    new_owner: Pubkey::new_from_array(new_owner),
                }
            }
            _ => return Err(ProgramError::InvalidInstructionData),
        })
    }
}

// Process create data ownership instruction
fn process_create_data_ownership(
    program_id: &Pubkey,
    accounts: &[AccountInfo],
    owner: Pubkey,
    data_hash: [u8; 32],
) -> ProgramResult {
    let account_info_iter = &mut accounts.iter();
    let data_account = next_account_info(account_info_iter)?;

    // Verify the account is owned by the program
    if data_account.owner != program_id {
        return Err(ProgramError::IncorrectProgramId);
    }

    // Initialize the data ownership account
    let mut data = data_account.data.borrow_mut();
    data[0..32].copy_from_slice(&owner.to_bytes());
    data[32..64].copy_from_slice(&data_hash);

    Ok(())
}

// Process transfer ownership instruction
fn process_transfer_ownership(
    program_id: &Pubkey,
    accounts: &[AccountInfo],
    new_owner: Pubkey,
) -> ProgramResult {
    let account_info_iter = &mut accounts.iter();
    let data_account = next_account_info(account_info_iter)?;
    let owner_account = next_account_info(account_info_iter)?;

    // Verify the account is owned by the program
    if data_account.owner != program_id {
        return Err(ProgramError::IncorrectProgramId);
    }

    // Verify the owner account is the signer
    if !owner_account.is_signer {
        return Err(ProgramError::MissingRequiredSignature);
    }

    // Update the ownership
    let mut data = data_account.data.borrow_mut();
    data[0..32].copy_from_slice(&new_owner.to_bytes());

    Ok(())
} 