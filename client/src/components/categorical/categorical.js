// jshint esversion: 6
import React from "react";
import _ from "lodash";
import {
  Button,
  Tooltip,
  Icon,
  ControlGroup,
  InputGroup
} from "@blueprintjs/core";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Category from "./category";

@connect(state => ({
  categoricalSelection: state.categoricalSelection
}))
class Categories extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      createAnnoModeActive: false
    };
  }

  handleCreateUserAnno = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "new user annotation category created"
    });
  };

  handleAnnoModeChange = () => {
    this.setState({ createAnnoModeActive: true });
  };

  render() {
    const { createAnnoModeActive } = this.state;
    const { categoricalSelection } = this.props;
    if (!categoricalSelection) return null;

    return (
      <div
        style={{
          padding: globals.leftSidebarSectionPadding
        }}
      >
        <p
          style={Object.assign({}, globals.leftSidebarSectionHeading, {
            marginTop: 4
          })}
        >
          Categorical Metadata
        </p>
        {_.map(
          categoricalSelection,
          (catState, catName) =>
            catName !== "free_annotation" ? (
              <Category
                key={catName}
                metadataField={catName}
                createAnnoModeActive={createAnnoModeActive}
              />
            ) : null
        )}
        <Category
          key={"free_annotation"}
          metadataField={"free_annotation"}
          isUserAnno
          createAnnoModeActive={createAnnoModeActive}
        />
        {!createAnnoModeActive ? (
          <Tooltip
            content="Create new categorical metadata field"
            position="bottom"
          >
            <Button
              data-testclass="createAnnoMode"
              data-testid="createAnnoMode"
              onClick={this.handleAnnoModeChange}
              intent="primary"
              icon="tag"
            >
              Create new categorical field
            </Button>
          </Tooltip>
        ) : (
          <ControlGroup>
            <InputGroup leftIcon="tag" />
            <Button
              data-testclass="createUserAnno"
              data-testid="createUserAnno"
              onClick={this.handleCreateUserAnno}
              intent="primary"
            >
              Add category
            </Button>
          </ControlGroup>
        )}
      </div>
    );
  }
}

export default Categories;
