use block_aligner::scan_block::*;
use block_aligner::scores::*;

use queries::{
    ScoresBAInsertInput
};

// cargo run --features simd_avx2 --release

fn main() {
    let query_result = run_query();
    let proteins_vec = &query_result.data.as_ref().unwrap().protein;

    // block-aligner
    // score
    let block_size = 256;
    let gaps = Gaps { open: -11, extend: -1 };
           
    // Iterate over the results of the cursor.
    for i in 0..proteins_vec.len() {
        let r = PaddedBytes::from_bytes::<AAMatrix>(&proteins_vec[i].sequence.as_bytes(), block_size);

        let mut scores = Vec::<ScoresBAInsertInput>::new();
        
        for j in 0..proteins_vec.len() {
            let q = PaddedBytes::from_bytes::<AAMatrix>(&proteins_vec[j].sequence.as_bytes(), block_size);

            // Align with traceback, but no x drop threshold.
            let a = Block::<_, true, false>::align(&q, &r, &BLOSUM62, gaps, block_size..=block_size, 0);
            let res = a.res();
            
            scores.push(ScoresBAInsertInput { compared_id: Some(proteins_vec[j].uniprot_id.clone()), reference_id: Some(proteins_vec[i].uniprot_id.clone()), score: Some(res.score), scope: Some("beta_casein".to_string())});
        }

        let mutation_result = run_mutation(scores.clone());
        println!("{:?}", mutation_result);

        println!("{}", i);
        println!("{}", proteins_vec[i].uniprot_id);

        // Hasura free tier 60 requests/minute
        use std::{thread, time};

        let one_second = time::Duration::from_millis(1000);
        let now = time::Instant::now();

        thread::sleep(one_second);

        assert!(now.elapsed() >= one_second);
    }
}

fn run_query() -> cynic::GraphQlResponse<queries::MyQuery> {
    use cynic::http::ReqwestBlockingExt;

    let query = build_query();

    reqwest::blocking::Client::new()
        .post("https://current-rabbit-87.hasura.app/v1/graphql")
        .run_graphql(query)
        .unwrap()
}

fn build_query() -> cynic::Operation<'static, queries::MyQuery> {
    use cynic::QueryBuilder;
    use queries::{
        MyQuery
    };

    MyQuery::build(())
}

fn run_mutation(scores: Vec<ScoresBAInsertInput>) -> cynic::GraphQlResponse<queries::MyMutation> {
    use cynic::http::ReqwestBlockingExt;

    let mutation = build_mutation(scores);

    reqwest::blocking::Client::new()
        .post("https://current-rabbit-87.hasura.app/v1/graphql")
        .run_graphql(mutation)
        .unwrap()
}

fn build_mutation(scores: Vec<ScoresBAInsertInput>) -> cynic::Operation<'static, queries::MyMutation> {
    use cynic::MutationBuilder;
    use queries::{
        MyMutation, MyMutationArguments
    };

    // vec![ScoresBAInsertInput { compared_id: Some("test".to_string()), reference_id: Some("test".to_string()), score: Some(10), scope: Some("beta_casein".to_string())}, ScoresBAInsertInput { compared_id: Some("test1".to_string()), reference_id: Some("test1".to_string()), score: Some(10), scope: Some("beta_casein".to_string())}]

    MyMutation::build(&MyMutationArguments {
        mutation_body: scores
    })
}

// https://generator.cynic-rs.dev
#[cynic::schema_for_derives(
    file = r#"schema.graphql"#,
    module = "schema",
)]
mod queries {
    use super::schema;

    #[derive(cynic::QueryFragment, Debug)]
    #[cynic(graphql_type = "query_root")]
    pub struct MyQuery {
        #[arguments(r#where = ProteinBoolExp { protein_type: Some(StringComparisonExp { _eq: Some("beta_casein".to_string()) }) })]
        pub protein: Vec<protein>,
    }

    #[derive(cynic::QueryFragment, Debug)]
    pub struct protein {
        pub uniprot_id: String,
        pub sequence: String,
    }

    #[derive(cynic::InputObject, Debug)]
    #[cynic(graphql_type = "protein_bool_exp")]
    pub struct ProteinBoolExp {
        #[cynic(rename = "protein_type")]
        pub protein_type: Option<StringComparisonExp>,
    }

    #[derive(cynic::InputObject, Debug)]
    #[cynic(graphql_type = "String_comparison_exp")]
    pub struct StringComparisonExp {
        #[cynic(rename = "_eq")]
        pub _eq: Option<String>,
    }

    // mutation
    #[derive(cynic::FragmentArguments, Debug)]
    pub struct MyMutationArguments {
        pub mutation_body: Vec<ScoresBAInsertInput>,
    }

    #[derive(cynic::QueryFragment, Debug)]
    #[cynic(
        graphql_type = "mutation_root",
        argument_struct = "MyMutationArguments"
    )]
    pub struct MyMutation {
        #[arguments(objects = args.mutation_body.clone())]
        pub insert_scores__ba: Option<scores_BA_mutation_response>,
    }

    #[derive(cynic::QueryFragment, Debug)]
    pub struct scores_BA_mutation_response {
        pub affected_rows: i32,
    }

    #[derive(cynic::InputObject, Clone, Debug)]
    #[cynic(graphql_type = "scores_BA_insert_input")]
    pub struct ScoresBAInsertInput {
        #[cynic(rename = "compared_id")]
        pub compared_id: Option<String>,
        #[cynic(rename = "reference_id")]
        pub reference_id: Option<String>,
        pub score: Option<i32>,
        pub scope: Option<String>
    }
}

mod schema {
    cynic::use_schema!("schema.graphql");
}

// cargo run --features simd_avx2 --release
