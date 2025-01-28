use musk::tracing::start_musk_tracing_subscriber;

fn main() {
    start_musk_tracing_subscriber();

    // debug!("This should be captured only by stdout");
    // info!("This should be captured only by stdout");
    // warn!("This should be captured only by stderr");
    // error!("This should be captured only by stderr");
}
