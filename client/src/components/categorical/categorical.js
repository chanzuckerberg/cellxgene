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
      createAnnoModeActive: false,
      newCategoryText: ""
    };
  }

  handleCreateUserAnno = () => {
    const { dispatch } = this.props;
    const { newCategoryText } = this.state;
    dispatch({
      type: "new user annotation category created",
      data: newCategoryText
    });
    this.setState({
      createAnnoModeActive: false,
      newCategoryText: ""
    });
  };

  handleAnnoModeChange = () => {
    console.log("anno mode changed");
    this.setState({ createAnnoModeActive: true });
  };

  handleCancelAnnoMode = () => {
    this.setState({ createAnnoModeActive: false });
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
              style={{ marginTop: 10 }}
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
          <form
            onSubmit={e => {
              e.preventDefault();
              this.handleCreateUserAnno();
            }}
          >
            <ControlGroup style={{ marginTop: 10 }}>
              <InputGroup
                autoFocus
                onChange={e =>
                  this.setState({ newCategoryText: e.target.value })
                }
                leftIcon="tag"
              />
              <Button
                data-testclass="createUserAnno"
                data-testid="createUserAnno"
                onClick={this.handleCreateUserAnno}
                intent="primary"
              >
                Add Category
              </Button>
              <Button
                data-testclass="cancelAnnoMode"
                data-testid="cancelAnnoMode"
                onClick={this.handleCancelAnnoMode}
                intent="default"
              >
                Cancel
              </Button>
            </ControlGroup>
          </form>
        )}
      </div>
    );
  }
}

export default Categories;
