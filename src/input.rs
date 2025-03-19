use std::io::{self, Write};

pub fn read_res() -> usize {
    loop {
        let mut input = String::new();
        print!("Enter resolution: ");
        io::stdout().flush().unwrap();
        io::stdin().read_line(&mut input).expect("Failed to read line.");

        match input.trim().parse::<usize>() {
            Ok(res) => return res,
            Err(_) => {
                println!("Please enter a valid integer");
            },
        }
    }
}

pub fn read_bounds() -> Result<String, Box<dyn std::error::Error>> {
    loop {
        println!("Enter boundary condition type (Dirichlet or Neumann):");
        io::stdout().flush()?;

        let mut input = String::new();
        io::stdin().read_line(&mut input)?;
        let mut input = input.trim().to_lowercase();

        if input == "d" {
            input.clear();
            input.push_str("dirichlet");
        } else if input == "n" {
            input.clear();
            input.push_str("neumann");
        }

        if input == "dirichlet" || input == "d" || input == "n" || input == "neumann" {
            println!("You entered: {}", input);
            return Ok(input);
        } else {
            println!("Invalid input. Please enter 'Dirichlet' or 'Neumann'.");
        }
    }
}
