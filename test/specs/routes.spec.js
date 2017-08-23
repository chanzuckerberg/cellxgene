import React from 'react';
import { shallow } from 'enzyme';
import { BrowserRouter as Router } from 'react-router-dom';

import Routes from '../../src/routes';

describe('Routes', () => {
  it('contains spec with an expectation', () => {
    expect(shallow(<Routes />).find(Router)).to.have.length(1);
  });
});
