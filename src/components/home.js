import React from 'react';
import { Link } from 'react-router-dom';
import _ from "lodash";
import Helmet from 'react-helmet';
import Container from './container';

import Categorical from "./categorical/categorical";
import Continuous from "./continuous/continuous";

const Home = () => {

  return (
    <Container>
      <Helmet title="cellxgene" />
      <h1><Link to="/page2">cellxgene</Link></h1>
      <Categorical/>
      <Continuous/>
    </Container>
  )
};

export default Home;

// <Categorical title={"Sample type"} category={types}/>
// <Categorical title={"Selection"} category={selection}/>
// <Categorical title={"Location"} category={location}/>
// <Categorical title={"Sample name"} category={names}/>
